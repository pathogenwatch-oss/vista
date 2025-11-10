import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import typer
from Bio.Seq import Seq

indel_finder = re.compile(r"\w-+\w")


def overlaps(coords1: tuple, coords2: tuple, threshold: int) -> bool:
    return min(coords1[1], coords2[1]) - max(coords1[0], coords2[0]) >= threshold


def find_frameshift(query: str, subject: str) -> bool:
    insert_matches = [
        insert for insert in indel_finder.findall(query) if len(insert) % 3 != 0
    ]
    if insert_matches:
        return True
    deletion_matches = [
        deletion for deletion in indel_finder.findall(subject) if len(deletion) % 3 != 0
    ]
    return bool(deletion_matches)


def find_premature_stop(dna: str, frame: int, includes_end: bool) -> bool:
    dna = dna.replace("-", "")[frame - 1 :]
    end_offset = (len(dna) % 3) * -1
    if end_offset != 0:
        includes_end = False
    dna = dna[0:end_offset]
    coding_seq = Seq(dna)
    translation = coding_seq.translate()
    if translation.count("*") < 1:
        return False
    else:
        end_check = 1 if includes_end else 0
        return str(translation).index("*") < len(translation) - end_check


def process_contig(
    contig_id: str, alignments: list[Any], lengths: dict[str, int], coverage: float
) -> dict[str, dict[str, Any]]:
    threshold = 60
    excluded = set()
    contig_keep = defaultdict(dict)

    for query_alignment in alignments:
        title = query_alignment.title.split(" ")[0]
        selected = []

        for hsp in query_alignment.hsps:
            name = f"{title}_{hsp.query_start}"
            if name in excluded:
                continue

            # Check for overlaps with all other HSPs
            for test_alignment in alignments:
                test_title = test_alignment.title.split(" ")[0]
                for test_hsp in test_alignment.hsps:
                    test_name = f"{test_title}_{test_hsp.query_start}"

                    if (
                        test_name in excluded
                        or name == test_name
                        or (hsp.sbjct_end - hsp.sbjct_start + 1) / lengths[title]
                        < coverage
                    ):
                        continue

                    if overlaps(
                        (hsp.query_start, hsp.query_end),
                        (test_hsp.query_start, test_hsp.query_end),
                        threshold,
                    ):
                        # Determine which HSP to exclude based on quality
                        hsp_score = hsp.align_length - hsp.gaps
                        test_score = test_hsp.align_length - test_hsp.gaps

                        if hsp_score < test_score:
                            excluded.add(test_name)
                        elif hsp.identities >= test_hsp.identities:
                            excluded.add(test_name)
                        else:
                            excluded.add(name)
                            break

            if name not in excluded:
                selected.append(hsp)

        if selected:
            contig_keep[title][contig_id] = selected

    return contig_keep


def select_matches(
    blast_records, lengths: dict[str, int], coverage: float
) -> dict[str, dict[str, Any]]:
    record_list = list(blast_records)
    kept: dict[str, dict[str, Any]] = dict()

    for contig_search in record_list:
        if len(contig_search.alignments) == 0:
            continue
        kept.update(
            process_contig(
                contig_search.query, contig_search.alignments, lengths, coverage
            )
        )
    return kept


def build_blastdb(
    db_dir: Path | str,
    name: str,
):
    fasta_path = db_dir / f"{name}.fasta.gz"
    db_path = db_dir / name
    # with gzip.open(fasta_path, "wb") as fasta_fh:
    #     for gene in genes:
    #         fasta_fh.write(sequences[gene].format("fasta").encode("utf-8"))

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
