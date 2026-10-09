from types import SimpleNamespace

from vista.blast_utils import (
    find_frameshift,
    find_premature_stop,
    overlaps,
    process_contig,
    select_matches,
)
from vista.files import read_sequences
from vista.search import classify_matches


def hsp(
    query_start: int,
    query_end: int,
    *,
    align_length: int,
    identities: int,
    gaps: int = 0,
    bits: float | None = None,
    sbjct_start: int = 1,
    sbjct_end: int = 100,
) -> SimpleNamespace:
    return SimpleNamespace(
        query_start=query_start,
        query_end=query_end,
        sbjct_start=sbjct_start,
        sbjct_end=sbjct_end,
        align_length=align_length,
        identities=identities,
        gaps=gaps,
        bits=align_length if bits is None else bits,
    )


def alignment(title: str, hit: SimpleNamespace) -> SimpleNamespace:
    return SimpleNamespace(title=title, hsps=[hit])


def test_overlapping_representatives_keep_higher_scoring_hit() -> None:
    lower_scoring = hsp(100, 199, align_length=80, identities=80)
    higher_scoring = hsp(120, 219, align_length=100, identities=90)

    kept = process_contig(
        "contig-1",
        [
            alignment("representative-low", lower_scoring),
            alignment("representative-high", higher_scoring),
        ],
        {"representative-low": 100, "representative-high": 100},
        coverage=0.8,
    )

    assert kept == {"representative-high": {"contig-1": [higher_scoring]}}


def test_overlapping_equal_score_representatives_keep_higher_identity_hit() -> None:
    lower_identity = hsp(100, 199, align_length=100, identities=80)
    higher_identity = hsp(120, 219, align_length=100, identities=90)

    kept = process_contig(
        "contig-1",
        [
            alignment("representative-low", lower_identity),
            alignment("representative-high", higher_identity),
        ],
        {"representative-low": 100, "representative-high": 100},
        coverage=0.8,
    )

    assert kept == {"representative-high": {"contig-1": [higher_identity]}}


def test_overlapping_representatives_prioritise_coverage_over_bitscore() -> None:
    lower_coverage = hsp(
        100,
        179,
        align_length=80,
        identities=80,
        bits=999,
        sbjct_end=80,
    )
    higher_coverage = hsp(120, 219, align_length=100, identities=90, bits=100)

    kept = process_contig(
        "contig-1",
        [
            alignment("representative-low", lower_coverage),
            alignment("representative-high", higher_coverage),
        ],
        {"representative-low": 100, "representative-high": 100},
        coverage=0.8,
    )

    assert kept == {"representative-high": {"contig-1": [higher_coverage]}}


def test_overlap_threshold_counts_inclusive_coordinates() -> None:
    assert overlaps((100, 199), (140, 239), threshold=60)
    assert not overlaps((100, 199), (141, 240), threshold=60)


def test_overlapping_representatives_prioritise_bitscore_over_identity() -> None:
    lower_bitscore = hsp(100, 199, align_length=100, identities=99, bits=100)
    higher_bitscore = hsp(120, 219, align_length=100, identities=90, bits=101)

    kept = process_contig(
        "contig-1",
        [
            alignment("representative-low", lower_bitscore),
            alignment("representative-high", higher_bitscore),
        ],
        {"representative-low": 100, "representative-high": 100},
        coverage=0.8,
    )

    assert kept == {"representative-high": {"contig-1": [higher_bitscore]}}


def test_overlapping_representatives_use_percent_identity_as_final_tiebreaker() -> None:
    lower_identity = hsp(100, 199, align_length=100, identities=90, bits=100)
    higher_identity = hsp(120, 219, align_length=100, identities=95, bits=100)

    kept = process_contig(
        "contig-1",
        [
            alignment("representative-low", lower_identity),
            alignment("representative-high", higher_identity),
        ],
        {"representative-low": 100, "representative-high": 100},
        coverage=0.8,
    )

    assert kept == {"representative-high": {"contig-1": [higher_identity]}}


def test_global_overlap_selection_retains_non_overlapping_flank_hit() -> None:
    middle = hsp(41, 140, align_length=100, identities=90, bits=100)
    left = hsp(1, 100, align_length=100, identities=90, bits=90)
    right = hsp(81, 180, align_length=100, identities=90, bits=110)

    kept = process_contig(
        "contig-1",
        [
            alignment("representative-middle", middle),
            alignment("representative-left", left),
            alignment("representative-right", right),
        ],
        {
            "representative-middle": 100,
            "representative-left": 100,
            "representative-right": 100,
        },
        coverage=0.8,
    )

    assert kept == {
        "representative-left": {"contig-1": [left]},
        "representative-right": {"contig-1": [right]},
    }


def test_partial_hit_is_discarded_even_when_it_is_the_only_alignment() -> None:
    partial_hit = hsp(100, 178, align_length=79, identities=79, sbjct_end=79)

    kept = process_contig(
        "contig-1",
        [alignment("representative", partial_hit)],
        {"representative": 100},
        coverage=0.8,
    )

    assert kept == {}


def test_reverse_strand_hit_uses_absolute_reference_span_for_coverage() -> None:
    reverse_hit = hsp(
        100,
        179,
        align_length=80,
        identities=80,
        sbjct_start=100,
        sbjct_end=21,
    )

    kept = process_contig(
        "contig-1",
        [alignment("representative", reverse_hit)],
        {"representative": 100},
        coverage=0.8,
    )

    assert kept == {"representative": {"contig-1": [reverse_hit]}}


def test_hsps_with_the_same_query_start_are_still_compared() -> None:
    lower_scoring = hsp(100, 199, align_length=80, identities=80)
    higher_scoring = hsp(100, 219, align_length=100, identities=90)

    kept = process_contig(
        "contig-1",
        [SimpleNamespace(title="representative", hsps=[lower_scoring, higher_scoring])],
        {"representative": 100},
        coverage=0.8,
    )

    assert kept == {"representative": {"contig-1": [higher_scoring]}}


def test_matches_from_multiple_contigs_are_retained() -> None:
    first_hit = hsp(1, 100, align_length=100, identities=100)
    second_hit = hsp(1, 100, align_length=100, identities=100)
    records = [
        SimpleNamespace(query="contig-1", alignments=[alignment("representative", first_hit)]),
        SimpleNamespace(query="contig-2", alignments=[alignment("representative", second_hit)]),
    ]

    kept = select_matches(records, {"representative": 100}, coverage=0.8)

    assert kept == {
        "representative": {
            "contig-1": [first_hit],
            "contig-2": [second_hit],
        }
    }


def test_frameshift_detection_uses_gap_length_only() -> None:
    assert find_frameshift("ATG-C", "ATGGC")
    assert not find_frameshift("ATG---C", "ATGAAAC")


def test_premature_stop_is_detected_in_full_codon_sequence() -> None:
    assert find_premature_stop("ATGTAAATG", frame=1, includes_end=False)
    assert not find_premature_stop("ATGAAATAA", frame=1, includes_end=True)


def test_complete_match_with_terminal_stop_is_not_disrupted() -> None:
    complete_hsp = SimpleNamespace(
        query_start=1,
        query_end=9,
        sbjct_start=1,
        sbjct_end=9,
        frame=(1, 1),
        query="ATGAAATAA",
        sbjct="ATGAAATAA",
        identities=9,
    )

    matches = classify_matches("representative", {"contig-1": [complete_hsp]}, 9)

    assert matches[0].isComplete
    assert not matches[0].isDisrupted


def test_partial_match_uses_reference_coding_phase_for_stop_detection() -> None:
    partial_hsp = SimpleNamespace(
        query_start=1,
        query_end=8,
        sbjct_start=2,
        sbjct_end=9,
        frame=(1, 1),
        query="TGAAATAA",
        sbjct="TGAAATAA",
        identities=8,
    )

    matches = classify_matches("representative", {"contig-1": [partial_hsp]}, 9)

    assert not matches[0].isDisrupted


def test_read_sequences_reads_an_uncompressed_fasta(tmp_path) -> None:
    (tmp_path / "references.fasta").write_text(">reference\nATG\n")

    sequences = read_sequences(tmp_path)

    assert str(sequences["reference"].seq) == "ATG"
