import io
import json
import multiprocessing
import os
import sys
import uuid
from multiprocessing.pool import ThreadPool as Pool
from typing import Any, Tuple

import typer
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

from vista.io import read_metadata, read_sequence_lengths, read_sequences, build_blastdb
from vista.search import serogroup_assignment, biotype_assignment, virulence_assignments, cluster_assignments
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
        out = os.path.join('/tmp', str(uuid.uuid4()))
        blastn_cline = NcbiblastnCommandline(query=query_fasta, db=blast_db, evalue=evalue, outfmt=5)
        stdout, stderr = blastn_cline()
        with io.StringIO(stdout) as f:
            selected_records = clean_matches(NCBIXML.parse(f), lengths, coverage)
        if name == 'virulenceGenes':
            return 'virulenceGenes', virulence_assignments(libraries[name], selected_records, lengths)
        elif name == 'biotypeMarkers':
            return 'biotypeMarkers', biotype_assignment(selected_records, lengths)
        elif name == 'serogroupMarkers':
            return 'serogroupMarkers', serogroup_assignment(libraries[name], selected_records, lengths)
        elif name == 'virulenceSets':
            return 'virulenceSets', cluster_assignments(libraries[name], selected_records, lengths)
        else:
            raise ValueError

    # Submit blasts in parallel
    vista_result = dict()
    with Pool(processes=cpus) as pool:
        for library, result in pool.map(library_search, list(libraries.keys())):
            vista_result[library] = result

    print(json.dumps(vista_result, default=lambda x: x.__dict__), file=sys.stdout)


@app.command()
def build(data_dir: str = 'data', db_dir: str = 'data'):
    metadata = read_metadata(data_dir)
    sequences = read_sequences(data_dir)
    libraries = metadata['libraries']

    # First build simple marker databases
    virulence_sets = 'virulenceSets'
    for name in libraries.keys():
        if name == virulence_sets:
            continue
        build_blastdb(data_dir, db_dir, [record["gene"] for record in libraries[name]['genes']], name, sequences)

    # Then build a DB for the virulence gene clusters
    virulence_genes = set()
    for cluster in libraries[virulence_sets]:
        virulence_genes.update(set(libraries[virulence_sets][cluster]["genes"]))
    build_blastdb(data_dir, db_dir, virulence_genes, virulence_sets, sequences)


if __name__ == '__main__':
    app()
