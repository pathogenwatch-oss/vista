import os
from typing import Dict

import toml
from Bio import SeqIO


def read_metadata(location: str) -> dict:
    metadata_path = os.path.join(location, 'metadata.toml')
    with open(metadata_path, 'r') as m_fh:
        return toml.load(m_fh)


def read_sequence_lengths(location: str) -> Dict[str, int]:
    fasta_path = os.path.join(location, 'references.fasta')
    with open(fasta_path, 'r') as f_fh:
        return {a: len(b) for a, b in SeqIO.to_dict(SeqIO.parse(f_fh, 'fasta', )).items()}