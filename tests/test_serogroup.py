from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from vista.vista import app
from vista.search import library_search


ROOT = Path(__file__).parent
RESOURCE_DIR = ROOT.parent / "src" / "vista" / "resources"
MARKERS_FASTA = RESOURCE_DIR / "serogroupMarkers.fasta.gz"
RUNNER = CliRunner()

SEROGROUP_LIBRARY = {
    "serogroupMarkers": {
        "genes": [
            {"name": "rfbV", "type": "O1"},
            {"name": "wbfZ", "type": "O139"},
        ]
    }
}


def marker_sequences() -> dict[str, str]:
    """Read the marker references used to make exact synthetic query FASTAs."""
    sequences: dict[str, str] = {}
    with gzip.open(MARKERS_FASTA, "rt") as fasta:
        name: str | None = None
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                name:str = line[1:]
                sequences[name] = ""
            elif name:
                sequences[name] += line
    return sequences


def search_serogroup(query_fasta: Path) -> dict:
    _, result = library_search(
        "serogroupMarkers",
        SEROGROUP_LIBRARY,
        str(RESOURCE_DIR),
        evalue=1e-20,
        coverage=0.8,
        query_fasta=query_fasta,
    )
    return result


def test_env_seawater_is_non_o1_o139() -> None:
    result = search_serogroup(ROOT / "resources" / "env_seawater.fasta")

    assert result == {
        "serogroup": "non-O1/O139",
        "serogroupMarkers": [
            {"name": "rfbV", "type": "O1", "matches": []},
            {"name": "wbfZ", "type": "O139", "matches": []},
        ],
    }


@pytest.mark.parametrize(
    ("markers", "expected_serogroup"),
    [
        (("rfbV",), "O1"),
        (("wbfZ",), "O139"),
        (("rfbV", "wbfZ"), "O1;O139"),
    ],
)
def test_serogroup_marker_combinations(
    tmp_path: Path, markers: tuple[str, ...], expected_serogroup: str
) -> None:
    sequences = marker_sequences()
    query_fasta = tmp_path / "synthetic_markers.fasta"
    query_fasta.write_text(
        "".join(f">{marker}\n{sequences[marker]}\n" for marker in markers)
    )

    result = search_serogroup(query_fasta)

    assert result["serogroup"] == expected_serogroup
    assert [marker["name"] for marker in result["serogroupMarkers"]] == [
        "rfbV",
        "wbfZ",
    ]
    assert [bool(marker["matches"]) for marker in result["serogroupMarkers"]] == [
        marker in markers for marker in ("rfbV", "wbfZ")
    ]


def without_version_fields(value):
    if isinstance(value, dict):
        return {
            key: without_version_fields(item)
            for key, item in value.items()
            if key != "version"
        }
    if isinstance(value, list):
        return [without_version_fields(item) for item in value]
    return value


def test_3247_cn_output_regression() -> None:
    result = RUNNER.invoke(
        app,
        ["search", str(ROOT / "resources" / "3247-CN.fasta"), "--cpus", "1"],
    )

    assert result.exit_code == 0, result.exception
    output = json.loads(result.stdout)
    normalized_output = json.dumps(
        without_version_fields(output), sort_keys=True
    ).encode("utf-8")

    actual_hash = hashlib.sha256(normalized_output).hexdigest()
    assert actual_hash == (
        "e7de562f2002e551a7ee035f660d6a0ecb4bf7a70e0a0bf49401e46e3bef7e06"
    ), actual_hash
