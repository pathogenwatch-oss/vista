import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

from Bio.Seq import Seq

indel_finder = re.compile(r"-+")


def overlaps(
    coords1: tuple[int, int], coords2: tuple[int, int], threshold: int
) -> bool:
    return min(coords1[1], coords2[1]) - max(coords1[0], coords2[0]) + 1 >= threshold


def find_frameshift(query: str, subject: str) -> bool:
    return any(
        len(indel) % 3 != 0
        for aligned_sequence in (query, subject)
        for indel in indel_finder.findall(aligned_sequence)
    )


def find_premature_stop(dna: str, frame: int, includes_end: bool) -> bool:
    dna = dna.replace("-", "")[frame - 1 :]
    remainder = len(dna) % 3
    if remainder:
        includes_end = False
        dna = dna[:-remainder]
    if not dna:
        return False
    coding_seq = Seq(dna)
    translation = coding_seq.translate()
    terminal_stop_offset = 1 if includes_end else 0
    return "*" in translation[: len(translation) - terminal_stop_offset]


def process_contig(
    contig_id: str, alignments: list[Any], lengths: dict[str, int], coverage: float
) -> dict[str, dict[str, Any]]:
    threshold = 60
    contig_keep = defaultdict(dict)

    def has_sufficient_coverage(alignment_hsp: Any, alignment_title: str) -> bool:
        return reference_coverage(alignment_hsp, alignment_title) >= coverage

    def reference_coverage(alignment_hsp: Any, alignment_title: str) -> float:
        return (
            abs(alignment_hsp.sbjct_end - alignment_hsp.sbjct_start) + 1
        ) / lengths[alignment_title]

    def quality(alignment_hsp: Any, alignment_title: str) -> tuple[float, float, float]:
        return (
            reference_coverage(alignment_hsp, alignment_title),
            alignment_hsp.bits,
            alignment_hsp.identities / alignment_hsp.align_length,
        )

    candidates = []
    for query_alignment in alignments:
        title = query_alignment.title.split(" ")[0]
        for hsp_index, hsp in enumerate(query_alignment.hsps):
            if not has_sufficient_coverage(hsp, title):
                continue
            candidates.append((title, hsp_index, hsp))

    selected = set()
    selected_hsps = []
    for title, hsp_index, candidate_hsp in sorted(
        candidates, key=lambda candidate: quality(candidate[2], candidate[0]), reverse=True
    ):
        if any(
            overlaps(
                (candidate_hsp.query_start, candidate_hsp.query_end),
                (selected_hsp.query_start, selected_hsp.query_end),
                threshold,
            )
            for selected_hsp in selected_hsps
        ):
            continue
        selected.add((title, hsp_index))
        selected_hsps.append(candidate_hsp)

    for title, hsp_index, candidate_hsp in candidates:
        if (title, hsp_index) in selected:
            contig_keep[title].setdefault(contig_id, []).append(candidate_hsp)

    return contig_keep


def select_matches(
    blast_records, lengths: dict[str, int], coverage: float
) -> dict[str, dict[str, Any]]:
    record_list = list(blast_records)
    kept: dict[str, dict[str, Any]] = dict()

    for contig_search in record_list:
        if len(contig_search.alignments) == 0:
            continue
        contig_matches = process_contig(
            contig_search.query, contig_search.alignments, lengths, coverage
        )
        for title, matches in contig_matches.items():
            kept.setdefault(title, {}).update(matches)
    return kept


def build_blastdb(
    db_dir: Path,
    name: str,
):
    fasta_path = db_dir / f"{name}.fasta.gz"
    db_path = db_dir / name

    gunzip_proc = subprocess.Popen(["gunzip", "-c", fasta_path], stdout=subprocess.PIPE)
    makeblastdb_proc = subprocess.Popen(
        [
            "makeblastdb",
            "-in",
            "-",  # Read from stdin
            "-title",
            name,
            "-out",
            str(db_path),
            "-dbtype",
            "nucl",
            "-parse_seqids",
        ],
        stdin=gunzip_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Allow gunzip to receive a SIGPIPE if makeblastdb exits.
    if gunzip_proc.stdout:
        gunzip_proc.stdout.close()

    stderr = makeblastdb_proc.communicate()[1]
    gunzip_proc.wait()

    if makeblastdb_proc.returncode != 0:
        print(f"Failed to build {name} database", file=sys.stderr)
        print(stderr.decode(), file=sys.stderr)
        raise ChildProcessError
