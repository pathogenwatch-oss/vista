import dataclasses
import hashlib
import io
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

from vista.blast_utils import find_frameshift, find_premature_stop, select_matches
from vista.files import read_sequence_lengths


@dataclasses.dataclass
class Match:
    queryId: str
    contigId: str
    queryStart: int
    queryEnd: int
    refStart: int
    refEnd: int
    frame: int
    isForward: bool
    isComplete: bool
    isDisrupted: bool
    isExact: bool
    identity: float

    # def __dict__(self) -> dict[str, Any]:
    #     """Custom dict representation excluding private fields"""
    #     return {
    #         k: v for k, v in dataclasses.asdict(self).items() if not k.startswith("_")
    #     }

    def overlaps(self, other: "Match", threshold: int) -> bool:
        """Check if this match overlaps with another match"""
        if self.contigId != other.contigId:
            return False

        # Calculate overlap between query coordinates
        overlap_start = max(self.queryStart, other.queryStart)
        overlap_end = min(self.queryEnd, other.queryEnd)
        overlap_length = max(0, overlap_end - overlap_start + 1)

        return overlap_length > threshold

    def extract(self, fasta: Dict[str, Any]) -> Seq:
        """Return the sequence of the match extracted from the FASTA file"""
        contig_seq = str(fasta[self.contigId].seq)

        # Extract sequence using 0-based indexing (BLAST coordinates are 1-based)
        start_idx = self.queryStart - 1
        end_idx = self.queryEnd
        extracted_seq = Seq(contig_seq[start_idx:end_idx])

        # Reverse complement if match is on reverse strand
        if not self.isForward:
            extracted_seq = extracted_seq.reverse_complement()

        return extracted_seq


@dataclasses.dataclass
class AaMatch(Match):
    aaIdentity: float
    aaName: str
    ntName: str

    @staticmethod
    def from_match(match: Match, aa_identity: float) -> "AaMatch":
        return AaMatch(
            **dataclasses.asdict(match),
            aaIdentity=aa_identity,
            aaName=match.queryId,
            ntName=match.queryId,
        )

    @staticmethod
    def calculate_aa_identity(
        dna_match: Match, query_seq: Seq, ref_seq: Seq
    ) -> "AaMatch":
        """Calculates the sequence identity between the query and reference sequences using Needleman-Wunsch global alignment"""
        # Translate sequences to amino acids
        query_aa = str(query_seq.translate())
        ref_aa = str(ref_seq.translate())

        if not query_aa or not ref_aa:
            aa_identity = 0.0
        else:
            # Create PairwiseAligner for Needleman-Wunsch global alignment
            aligner = PairwiseAligner()
            aligner.mode = "global"
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5

            # Perform alignment
            alignments = aligner.align(query_aa, ref_aa)

            if not alignments:
                aa_identity = 0.0
            else:
                # Take the best alignment
                best_alignment = alignments[0]
                # Count identical positions (excluding gaps)
                identical_positions = 0
                for q, r in zip(best_alignment[0], best_alignment[1]):
                    if q == r and q != "-" and r != "-":
                        identical_positions += 1

                # Calculate identity as percentage of the longer original sequence
                max_len = max(len(query_aa), len(ref_aa))
                aa_identity = round((identical_positions / max_len) * 100, 2)

        if aa_identity == 100.0:
            aa_name = dna_match.queryId
        else:
            aa_name = hashlib.sha1(str(ref_aa).encode("utf-8")).hexdigest()

        # Create AaMatch with all fields from the original Match plus aa_identity
        return AaMatch(
            **dataclasses.asdict(dna_match),
            aaIdentity=aa_identity,
            aaName=aa_name,
            ntName=hashlib.sha1(str(query_seq).encode("utf-8")).hexdigest(),
        )


def classify_matches(
    query_id: str, family_matches: Dict[str, Any], ref_length: int
) -> List[Match]:
    matches = list()
    for contig_id, hsps in family_matches.items():
        for hsp in hsps:
            ref_start, ref_end = sorted([hsp.sbjct_start, hsp.sbjct_end])
            is_forward = 0 < hsp.frame[1]
            is_complete = ref_start == 1 and ref_end == ref_length
            is_exact = is_complete and hsp.query == hsp.sbjct
            query_seq = (
                hsp.query if is_forward else str(Seq(hsp.query).reverse_complement())
            )
            is_disrupted = find_frameshift(query_seq, hsp.sbjct) or find_premature_stop(
                query_seq, hsp.frame[0], ref_end < ref_length
            )
            percent_identity = round(
                (float(hsp.identities) / float(ref_length)) * 100, 2
            )
            matches.append(
                Match(
                    query_id,
                    contig_id,
                    hsp.query_start,
                    hsp.query_end,
                    ref_start,
                    ref_end,
                    hsp.frame[0],
                    is_forward,
                    is_complete,
                    is_disrupted,
                    is_exact,
                    percent_identity,
                )
            )
    return matches


def determine_status(matches: List[Match]) -> str:
    if not matches:
        return "Not found"

    return (
        "Present"
        if any(not match.isDisrupted and match.isComplete for match in matches)
        else "Incomplete"
    )


def gather_hits_by_name(
    records: Dict[str, Any], marker_name: str, lengths: Dict[str, int]
) -> List[Match]:
    return [
        hit
        for hit in classify_matches(
            marker_name, records[marker_name], lengths[marker_name]
        )
        if not hit.isDisrupted or hit.isComplete
    ]


def at_least_one_exact_match(hits: List[Match]) -> bool:
    return 0 < len(only_exact_matches(hits))


def only_exact_matches(hits: List[Match]) -> List[Match]:
    return [hit for hit in hits if hit.isExact]


def serogroup_assignment(
    library: Dict[str, Any], records, lengths: Dict[str, int]
) -> Dict[str, Any]:
    type_markers = list()
    types = list()
    for marker in library["genes"]:
        if marker["name"] in records.keys():
            hits = gather_hits_by_name(records, marker["name"], lengths)

            if at_least_one_exact_match(hits):
                types.append(marker["type"])
            marker["matches"] = hits
        else:
            marker["matches"] = []
        type_markers.append(marker)

    tag = ";".join(types)
    return {"serogroup": tag, "serogroupMarkers": type_markers}


def virulence_assignments(
    library: Dict[str, Any], records: Dict[str, Any], lengths: Dict[str, int]
) -> Dict[str, Any]:
    result = list()

    for gene in library["genes"]:
        if gene["name"] in records.keys():
            virulence_hits = classify_matches(
                gene["name"], records[gene["name"]], lengths[gene["name"]]
            )
            gene["status"] = determine_status(virulence_hits)
            gene["matches"] = virulence_hits
        else:
            gene["status"] = "Not found"
            gene["matches"] = []
        result.append(gene)
    return {"virulenceGenes": result}


def cluster_assignments(
    library: Dict[str, Any], records: Dict[str, Any], lengths: Dict[str, int]
) -> Dict[str, Any]:
    result = list()

    for cluster_id, cluster in library.items():
        cluster["id"] = cluster_id
        cluster["matches"] = dict()
        for gene in cluster["genes"]:
            gene_profile = dict()
            if gene in records.keys():
                cluster_hits = classify_matches(gene, records[gene], lengths[gene])
                gene_profile["status"] = determine_status(cluster_hits)
                gene_profile["matches"] = cluster_hits
            else:
                gene_profile["status"] = "Not found"
                gene_profile["matches"] = []
            cluster["matches"][gene] = gene_profile
        cluster["present"] = [
            gene
            for gene, profile in cluster["matches"].items()
            if "Present" == profile["status"]
        ]
        cluster["missing"] = [
            gene
            for gene, profile in cluster["matches"].items()
            if "Not found" == profile["status"]
        ]
        cluster["incomplete"] = [
            gene
            for gene, profile in cluster["matches"].items()
            if "Incomplete" == profile["status"]
        ]
        # cluster['complete'] = len(cluster['present']) == len(cluster['genes'])
        cluster["status"] = (
            "Present"
            if len(cluster["present"]) == len(cluster["genes"])
            else "Incomplete"
            if len(cluster["present"]) > 0
            else "Not found"
        )
        result.append(cluster)
    return {"virulenceClusters": result}


def library_search(
    name: str,
    libraries: dict[str, Any],
    data_path: str,
    evalue: float,
    coverage: float,
    query_fasta: Path | str,
    num_threads: int = 1,
) -> Tuple[str, Any]:
    print(f"Searching {name} library", file=sys.stderr)
    blast_db = os.path.join(data_path, name)
    lengths = read_sequence_lengths(os.path.join(data_path, f"{name}.fasta.gz"))
    # Construct the blastn command
    # Run the blastn command
    blast_res = subprocess.run(
        [
            "blastn",
            "-query",
            str(query_fasta),
            "-db",
            blast_db,
            "-evalue",
            str(evalue),
            "-outfmt",
            "5",  # XML output format
            "-num_threads",
            str(num_threads),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    with io.StringIO(blast_res.stdout) as f:
        selected_records = select_matches(NCBIXML.parse(f), lengths, coverage)
    if name == "virulenceGenes":
        return "virulenceGenes", virulence_assignments(
            libraries[name], selected_records, lengths
        )
    elif name == "serogroupMarkers":
        return "serogroupMarkers", serogroup_assignment(
            libraries[name], selected_records, lengths
        )
    elif name == "virulenceClusters":
        return "virulenceClusters", cluster_assignments(
            libraries[name], selected_records, lengths
        )
    else:
        print(f"The BLAST search for {name} library was unsuccessful.")
        raise ValueError
