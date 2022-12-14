import json
import os
import sys
import uuid
from typing import List, Dict

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq

from vista.io import read_metadata, read_sequence_lengths
from vista.vista_utils import find_frameshift, find_premature_stop, remove_overlaps


class Match:
    def __init__(self, contig_id: str, query_start: int, query_end: int, ref_start: int, ref_end: int, frame: int,
                 is_forward: bool, is_complete: bool, is_disrupted: bool, is_exact: bool, identity: float):
        self.contigId = contig_id
        self.queryStart = query_start
        self.queryEnd = query_end
        self.refStart = ref_start
        self.refEnd = ref_end
        self.frame = frame
        self.isForward = is_forward
        self.isComplete = is_complete
        self.isDisrupted = is_disrupted
        self.isExact = is_exact
        self.identity = identity


def classify_matches(family_matches: dict, ref_length: int) -> List[Match]:
    matches = list()
    for contig_id, hsps in family_matches.items():
        for hsp in hsps:
            ref_start, ref_end = sorted([hsp.sbjct_start, hsp.sbjct_end])
            is_forward = 0 < hsp.frame[1]
            is_complete = ref_start == 1 and ref_end == ref_length
            is_exact = is_complete and hsp.query == hsp.sbjct
            query_seq = hsp.query if is_forward else str(Seq(hsp.query).reverse_complement())
            is_disrupted = find_frameshift(query_seq, hsp.sbjct) or find_premature_stop(query_seq, hsp.frame[0],
                                                                                        ref_end < ref_length)
            percent_identity = round((float(hsp.identities) / float(ref_length)) * 100, 2)
            matches.append(
                Match(contig_id, hsp.query_start, hsp.query_end, ref_start, ref_end, hsp.frame[0],
                      is_forward, is_complete, is_disrupted, is_exact, percent_identity))
    return matches


# Considered present if any match is complete and not disrupted.
def determine_status(matches: List[Match]) -> str:
    if len(matches) == 0:
        return 'Not found'
    status = 'Present'
    for match in matches:
        if match.isDisrupted or not match.isComplete:
            status = 'Incomplete'
        else:
            return 'Present'
    return status


def gather_hits_for_gene(records, marker_name):
    return [hit for hit in classify_matches(records[marker_name], lengths[marker_name]) if
            not hit.isDisrupted or not hit.isComplete]


def at_least_one_exact_match(hits):
    return 0 < len(only_exact_matches(hits))


def only_exact_matches(hits):
    return [hit for hit in hits if hit.isExact]


def extract_type(markers: List, records) -> (str, List):
    type_markers = list()
    types = list()
    for marker in markers:
        if marker['gene'] in records.keys():
            hits = gather_hits_for_gene(records, marker['gene'])

            if at_least_one_exact_match(hits):
                types.append(marker['name'])
            # else:
            #     types.append(marker['name'] + '*')
            marker['matches'] = hits
        else:
            marker['matches'] = []
        type_markers.append(marker)

    tag = ';'.join(types)
    return tag, type_markers


def extract_biotype(records: Dict) -> (str, List):
    types = list()
    type_markers = list()

    if 'wbfZ' in records.keys():
        hits = only_exact_matches(gather_hits_for_gene(records, 'wbfZ'))
        if 0 < len(hits):
            types.append('O139')
            type_markers.append({'name': 'O139', 'gene': 'wbfZ', 'matches': hits})

    rfbV_present = 'rfbV' in records.keys()
    if rfbV_present:
        rfbv_hits = only_exact_matches(gather_hits_for_gene(records, 'rfbV'))
        if 0 < len(rfbv_hits):
            type_markers.append({'name': 'O1', 'gene': 'rfbV', 'matches': rfbv_hits})
            modern_o1 = False
            has_ctxb = False
            all_ctxb_matches = list()
            modern_o1_schema = {'ctxB1': 'O1 classical', 'ctxB3': 'O1 El Tor', 'ctxB7': 'O1 Haiti'}
            for allele in modern_o1_schema:
                if allele in records.keys():
                    hits = gather_hits_for_gene(records, allele)
                    has_ctxb = True
                    all_ctxb_matches.extend(hits)
                    if at_least_one_exact_match(hits):
                        modern_o1 = True
                        types.append(modern_o1_schema[allele])
                        type_markers.append({'name': modern_o1_schema[allele], 'gene': allele, 'matches': hits})
            if not modern_o1 and has_ctxb:
                types.append('O1 pathogenic')
                type_markers.append({'name': 'O1 pathogenic', 'gene': 'ctxB', 'matches': all_ctxb_matches})
            if not modern_o1 and not has_ctxb:
                types.append('O1 environmental')
                type_markers.append({'name': 'O1 pathogenic', 'gene': 'ctxB', 'matches': all_ctxb_matches})

    return ';'.join(types), type_markers


query_fasta = sys.argv[1]
data_path = sys.argv[2] if len(sys.argv) == 3 else sys.path[0] + '/data/'
out = str(uuid.uuid4())

blast_db = os.path.join(data_path, 'references')
metadata = read_metadata(data_path)
lengths = read_sequence_lengths(data_path)

# Run Blast against input FASTA
evalue = '1e-20'
cutoff = '0.8'

blastn_cline = NcbiblastnCommandline(query=query_fasta, db=blast_db, evalue=evalue, outfmt=5, out=out)
stdout, stderr = blastn_cline()
with open(out, 'r') as r_fh:
    selected_records = remove_overlaps(NCBIXML.parse(r_fh))

os.remove(out)
# Check library against matches
# Simple virulence makers

result = dict()
result['virulenceGenes'] = list()

for gene in metadata['virulence_genes']:
    if gene['name'] in selected_records.keys():
        virulence_hits = classify_matches(selected_records[gene['name']], lengths[gene['name']])
        gene['status'] = determine_status(virulence_hits)
        gene['matches'] = virulence_hits
    elif 'references' in gene.keys() and len(selected_records.keys() & set(gene['references'])) > 0:
        # will only ever be one reference matched
        reference = list(selected_records.keys() & set(gene['references']))[0]
        virulence_hits = classify_matches(selected_records[reference], lengths[reference])
        gene['status'] = determine_status(virulence_hits)
        gene['matches'] = virulence_hits
    else:
        gene['status'] = 'Not found'
        gene['matches'] = []
    result['virulenceGenes'].append(gene)

# Virulence clusters/operons
result['virulenceClusters'] = list()

for cluster_id, cluster in metadata['virulence_sets'].items():
    cluster['id'] = cluster_id
    cluster['matches'] = dict()
    for gene in cluster['genes']:
        gene_profile = dict()
        if gene in selected_records.keys():
            cluster_hits = classify_matches(selected_records[gene], lengths[gene])
            gene_profile['status'] = determine_status(cluster_hits)
            gene_profile['matches'] = cluster_hits
        else:
            gene_profile['status'] = 'Not found'
            gene_profile['matches'] = []
        cluster['matches'][gene] = gene_profile
    cluster['present'] = [gene for gene, profile in cluster['matches'].items() if 'Present' == profile['status']]
    cluster['missing'] = [gene for gene, profile in cluster['matches'].items() if 'Not found' == profile['status']]
    cluster['incomplete'] = [gene for gene, profile in cluster['matches'].items() if 'Incomplete' == profile['status']]
    # cluster['complete'] = len(cluster['present']) == len(cluster['genes'])
    cluster['status'] = 'Present' if len(cluster['present']) == len(cluster['genes']) else 'Incomplete' if len(
        cluster['present']) > 0 else 'Not found'
    result['virulenceClusters'].append(cluster)

# Serogroups & biotypes
result['serogroup'], result['serogroupMarkers'] = extract_type(metadata['serogroup_markers'], selected_records)
result['biotype'], result['biotypeMarkers'] = extract_biotype(selected_records)

print(json.dumps(result, default=lambda x: x.__dict__), file=sys.stdout)
