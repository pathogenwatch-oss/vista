from typing import List, Dict, Any

from Bio.Seq import Seq

from vista.utils import find_frameshift, find_premature_stop


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


def classify_matches(family_matches: Dict[str, Any], ref_length: int) -> List[Match]:
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


def gather_hits_for_gene(records: Dict[str, Any], marker_name: str, lengths: Dict[str, int]) -> List[Match]:
    return [hit for hit in classify_matches(records[marker_name], lengths[marker_name]) if
            not hit.isDisrupted or not hit.isComplete]


def at_least_one_exact_match(hits: List[Match]) -> bool:
    return 0 < len(only_exact_matches(hits))


def only_exact_matches(hits: List[Match]) -> List[Match]:
    return [hit for hit in hits if hit.isExact]


def serogroup_assignment(library: Dict[str, Any], records, lengths: Dict[str, int]) -> (str, List[Dict[str, str]]):
    type_markers = list()
    types = list()
    for marker in library['genes']:
        if marker['name'] in records.keys():
            hits = gather_hits_for_gene(records, marker['name'], lengths)

            if at_least_one_exact_match(hits):
                types.append(marker['type'])
            # else:
            #     types.append(marker['name'] + '*')
            marker['matches'] = hits
        else:
            marker['matches'] = []
        type_markers.append(marker)

    tag = ';'.join(types)
    return {"serogroup": tag, "serogroupMarkers": type_markers}


def biotype_assignment(records: Dict[str, Any], lengths: Dict[str, int]) -> (str, List[Dict[str, Any]]):
    types = list()
    type_markers = list()

    if 'wbfZ' in records.keys():
        hits = only_exact_matches(gather_hits_for_gene(records, 'wbfZ', lengths))
        if 0 < len(hits):
            types.append('O139')
            type_markers.append({'type': 'O139', 'name': 'wbfZ', 'matches': hits})

    rfbv_present = 'rfbV' in records.keys()
    if rfbv_present:
        rfbv_hits = only_exact_matches(gather_hits_for_gene(records, 'rfbV', lengths))
        if 0 < len(rfbv_hits):
            type_markers.append({'type': 'O1', 'name': 'rfbV', 'matches': rfbv_hits})
            modern_o1 = False
            has_ctxb = False
            all_ctxb_matches = list()
            modern_o1_schema = {'ctxB1': 'O1 classical', 'ctxB3': 'O1 El Tor', 'ctxB7': 'O1 Haiti'}
            for allele in modern_o1_schema:
                if allele in records.keys():
                    hits = gather_hits_for_gene(records, allele, lengths)
                    has_ctxb = True
                    all_ctxb_matches.extend(hits)
                    if at_least_one_exact_match(hits):
                        modern_o1 = True
                        types.append(modern_o1_schema[allele])
                        type_markers.append({'type': modern_o1_schema[allele], 'name': allele, 'matches': hits})
            if not modern_o1 and has_ctxb:
                types.append('O1 pathogenic')
                type_markers.append({'type': 'O1 pathogenic', 'name': 'ctxB', 'matches': all_ctxb_matches})
            if not modern_o1 and not has_ctxb:
                types.append('O1 environmental')
                type_markers.append({'type': 'O1 pathogenic', 'name': 'ctxB', 'matches': all_ctxb_matches})

    return {"biotype": ';'.join(types), "biotypeMarkers": type_markers}


def virulence_assignments(library: Dict[str, Any], records: Dict[str, Any], lengths: Dict[str, int]) -> (
        str, List[Dict[str, Any]]):
    result = list()

    for gene in library['genes']:
        if gene['name'] in records.keys():
            virulence_hits = classify_matches(records[gene['name']], lengths[gene['name']])
            gene['status'] = determine_status(virulence_hits)
            gene['matches'] = virulence_hits
        else:
            gene['status'] = 'Not found'
            gene['matches'] = []
        result.append(gene)
    return {"virulenceGenes": result}


def cluster_assignments(library: Dict[str, Any], records: Dict[str, Any], lengths: Dict[str, int]) -> (
        str, List[Dict[str, Any]]):
    result = list()

    for cluster_id, cluster in library.items():
        cluster['id'] = cluster_id
        cluster['matches'] = dict()
        for gene in cluster['genes']:
            gene_profile = dict()
            if gene in records.keys():
                cluster_hits = classify_matches(records[gene], lengths[gene])
                gene_profile['status'] = determine_status(cluster_hits)
                gene_profile['matches'] = cluster_hits
            else:
                gene_profile['status'] = 'Not found'
                gene_profile['matches'] = []
            cluster['matches'][gene] = gene_profile
        cluster['present'] = [gene for gene, profile in cluster['matches'].items() if 'Present' == profile['status']]
        cluster['missing'] = [gene for gene, profile in cluster['matches'].items() if 'Not found' == profile['status']]
        cluster['incomplete'] = [gene for gene, profile in cluster['matches'].items() if
                                 'Incomplete' == profile['status']]
        # cluster['complete'] = len(cluster['present']) == len(cluster['genes'])
        cluster['status'] = 'Present' if len(cluster['present']) == len(cluster['genes']) else 'Incomplete' if len(
            cluster['present']) > 0 else 'Not found'
        result.append(cluster)
    return {"virulenceClusters": result}
