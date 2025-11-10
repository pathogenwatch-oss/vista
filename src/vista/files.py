import gzip
from pathlib import Path
from typing import Any, Dict

import toml
from Bio import SeqIO


def read_metadata(location: Path) -> Dict[str, Any]:
    with open(location, "r") as m_fh:
        return toml.load(m_fh)


def read_sequence_lengths(fasta_path: Path|str) -> Dict[str, int]:
    with gzip.open(fasta_path, "rt") as f_fh:
        return {
            a: len(b.seq)
            for a, b in SeqIO.to_dict(
                SeqIO.parse(
                    f_fh,
                    "fasta",
                )
            ).items()
        }


def read_sequences(location: Path) -> Dict[str, SeqIO.SeqRecord]:
    fasta_path = location / "references.fasta"
    with gzip.open(fasta_path, "rb") as f_fh:
        return SeqIO.to_dict(
            SeqIO.parse(
                f_fh,
                "fasta",
            )
        )
