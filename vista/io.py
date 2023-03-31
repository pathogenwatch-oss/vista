import os
import subprocess
from typing import Dict, Any, Iterable

import toml
import typer
from Bio import SeqIO


def read_metadata(location: str) -> Dict[str, Any]:
    metadata_path = os.path.join(location, 'metadata.toml')
    with open(metadata_path, 'r') as m_fh:
        return toml.load(m_fh)


def read_sequence_lengths(fasta_path: str) -> Dict[str, int]:
    with open(fasta_path, 'r') as f_fh:
        return {a: len(b) for a, b in SeqIO.to_dict(SeqIO.parse(f_fh, 'fasta', )).items()}


def read_sequences(location: str) -> Dict[str, Any]:
    fasta_path = os.path.join(location, 'references.fasta')
    with open(fasta_path, 'r') as f_fh:
        return SeqIO.to_dict(SeqIO.parse(f_fh, 'fasta', ))


def build_blastdb(data_dir: str, db_dir: str, genes: Iterable[str], name: str, sequences: Dict[str, Any]):
    fasta_path = os.path.join(data_dir, f'{name}.fasta')
    db_path = os.path.join(db_dir, name)
    with open(fasta_path, 'w') as fasta_fh:
        for gene in genes:
            print(sequences[gene].format("fasta"), file=fasta_fh)
    process = subprocess.run(["makeblastdb", "-in", fasta_path, "-out", db_path, "-dbtype", "nucl", "-parse_seqids"])
    if process.returncode != 0:
        typer.echo(f'Failed to build {name} database')
        raise ChildProcessError
