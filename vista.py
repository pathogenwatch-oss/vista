import io
import json
import multiprocessing
import os
import subprocess
import sys
from multiprocessing.pool import ThreadPool as Pool
from typing import Any, Tuple

import typer
from Bio.Blast import NCBIXML

from vista.files import read_metadata, read_sequence_lengths, read_sequences, build_blastdb
from vista.search import serogroup_assignment, virulence_assignments, cluster_assignments
from vista.utils import clean_matches

app = typer.Typer()


@app.command()
def search(query_fasta: str, data_path: str = 'data', cpus: int = multiprocessing.cpu_count()):
    metadata = read_metadata(data_path)
    blast_defaults = metadata['defaults']['blast']
    libraries = metadata['libraries']
    evalue = blast_defaults['evalue']
    coverage = blast_defaults['coverage']

    def library_search(name: str) -> Tuple[str, Any]:
        blast_db = os.path.join(data_path, name)
        lengths = read_sequence_lengths(os.path.join(data_path, f'{name}.fasta'))
        # Construct the blastn command
        # Run the blastn command
        blast_res = subprocess.run(
            [
                "blastn",
                "-query",
                query_fasta,
                "-db",
                blast_db,
                "-evalue",
                str(evalue),
                "-outfmt",
                "5",  # XML output format
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        with io.StringIO(blast_res.stdout) as f:
            selected_records = clean_matches(NCBIXML.parse(f), lengths, coverage)
        if name == 'virulenceGenes':
            return 'virulenceGenes', virulence_assignments(libraries[name], selected_records, lengths)
        elif name == 'serogroupMarkers':
            return 'serogroupMarkers', serogroup_assignment(libraries[name], selected_records, lengths)
        elif name == 'virulenceClusters':
            return 'virulenceClusters', cluster_assignments(libraries[name], selected_records, lengths)
        else:
            raise ValueError

    # Submit blasts in parallel
    vista_result = dict()
    with Pool(processes=cpus) as pool:
        for library, result in pool.map(library_search, list(libraries.keys())):
            vista_result = vista_result | result

    print(json.dumps(vista_result, default=lambda x: x.__dict__), file=sys.stdout)


@app.command()
def build(data_dir: str = 'data', db_dir: str = 'data'):
    metadata = read_metadata(data_dir)
    sequences = read_sequences(data_dir)
    libraries = metadata['libraries']

    # First build simple marker databases
    virulence_sets = 'virulenceClusters'
    for name in libraries.keys():
        if name == virulence_sets:
            continue
        build_blastdb(data_dir, db_dir, [record["name"] for record in libraries[name]['genes']], name, sequences)

    # Then build a DB for the virulence gene clusters
    virulence_genes = set()
    for cluster in libraries[virulence_sets]:
        virulence_genes.update(set(libraries[virulence_sets][cluster]["genes"]))
    build_blastdb(data_dir, db_dir, virulence_genes, virulence_sets, sequences)


if __name__ == '__main__':
    app()
