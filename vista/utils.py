import re
from collections import defaultdict
from typing import Any, Dict, List

from Bio.Seq import Seq

indel_finder = re.compile(r'\w-+\w')


def overlaps(coords1: tuple, coords2: tuple, threshold: int):
    return min(coords1[1], coords2[1]) - max(coords1[0], coords2[0]) >= threshold


def find_frameshift(query: str, sbjct: str) -> bool:
    insert_matches = [insert for insert in indel_finder.findall(query) if len(insert) % 3 != 0]
    if len(insert_matches) != 0:
        return True
    deletion_matches = [deletion for deletion in indel_finder.findall(sbjct) if len(deletion) % 3 != 0]
    if len(deletion_matches) != 0:
        return True


def find_premature_stop(dna: str, frame: int, includes_end: bool) -> bool:
    dna = dna.replace('-', '')[frame - 1:]
    end_offset = (len(dna) % 3) * -1
    if end_offset != 0:
        includes_end = False
    dna = dna[0:end_offset]
    coding_seq = Seq(dna)
    translation = coding_seq.translate()
    if translation.count('*') < 1:
        return False
    else:
        end_check = 1 if includes_end else 0
        return str(translation).index('*') < len(translation) - end_check


def process_contig(contig_id: str, alignments: List[Any], lengths: Dict[str, int],
                   coverage: float) -> Dict[str, Dict[str, Any]]:
    threshold = 60

    excluded = set()
    contig_keep = defaultdict(dict)

    for query in range(0, len(alignments)):
        selected = list()
        query_alignment = alignments[query]

        title = query_alignment.title.split(' ')[0]

        for hsp in query_alignment.hsps:
            name = title + '_' + str(hsp.query_start)
            if name in excluded:
                continue
            for test in range(query, len(alignments)):
                test_ali = alignments[test]
                test_title = test_ali.title.split(' ')[0]
                for test_hsp in test_ali.hsps:
                    test_name = test_title + '_' + str(test_hsp.query_start)
                    if test_name in excluded:
                        continue
                    if name == test_name:
                        continue
                    if (hsp.sbjct_end - hsp.sbjct_start + 1) / lengths[title] < coverage:
                        continue
                    if overlaps((hsp.query_start, hsp.query_end), (test_hsp.query_start, test_hsp.query_end),
                                threshold):
                        if hsp.align_length - hsp.gaps < test_hsp.align_length - test_hsp.gaps:
                            excluded.add(test_name)
                            continue
                        elif hsp.identities >= test_hsp.identities:
                            excluded.add(test_name)
                        else:
                            excluded.add(name)
                            break
                    else:
                        continue
            if name not in excluded:
                selected.append(hsp)
        if len(selected) != 0:
            contig_keep[title][contig_id] = selected
    return contig_keep


def clean_matches(blast_records, lengths: Dict[str, int], coverage: float) -> Dict[str, Dict[str, Any]]:
    record_list = list(blast_records)
    kept = dict()

    for contig_search in record_list:
        if len(contig_search.alignments) == 0:
            continue
        kept.update(process_contig(contig_search.query, contig_search.alignments, lengths, coverage))
    return kept
