import json
import multiprocessing
import sys
from functools import partial
from importlib import resources
from multiprocessing.pool import ThreadPool as Pool
from pathlib import Path
from typing import Annotated, Any

import typer

from vista.blast_utils import build_blastdb
from vista.files import (
    read_metadata,
)
from vista.search import library_search

VIRULENCE_SETS = "virulenceClusters"

app = typer.Typer()


def search_parallelism(cpus: int, library_count: int) -> tuple[int, int]:
    workers = min(cpus, library_count)
    return workers, max(1, cpus // workers)


def package_resource_path(resource_name: str) -> Path:
    resource = resources.files("vista").joinpath(resource_name)
    if isinstance(resource, Path):
        return resource
    raise RuntimeError(
        "VISTA resources must be installed as filesystem files; "
        "archive-based package resources cannot be used for BLAST databases."
    )


@app.command()
def search(
    query_fasta: Annotated[
        Path,
        typer.Argument(
            help="The path to the query FASTA file.",
            file_okay=True,
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    metadata_toml: Annotated[
        Path | None,
        typer.Option(
            "-m",
            "--metadata-toml",
            help="The path to the metadata TOML file. Defaults to './src/vista/config/metadata.toml' or package resources.",
            file_okay=True,
            dir_okay=False,
            readable=True,
        ),
    ] = None,
    resources_path: Annotated[
        Path | None,
        typer.Option(
            "-d",
            "--data-path",
            help="The location of the BLAST databases. Defaults to './src/vista/resources' or package resources.",
            file_okay=False,
            dir_okay=True,
            readable=True,
        ),
    ] = None,
    cpus: Annotated[
        int,
        typer.Option(
            "-c",
            "--cpus",
            help="The number of processes to run simultaneously. Defaults to the number of CPUs.",
        ),
    ] = multiprocessing.cpu_count(),
):
    if metadata_toml is None:
        if Path("pyproject.toml").exists():
            metadata_toml = Path("src/vista/config/metadata.toml")
        else:
            metadata_toml = package_resource_path("config/metadata.toml")

    if resources_path is None:
        resources_path = get_resources_dir()

    if not metadata_toml.is_file():
        print(f"Error: metadata_toml '{metadata_toml}' does not exist.", file=sys.stderr)
        raise typer.Exit(code=1)

    if not resources_path.is_dir():
        print(f"Error: data_path '{resources_path}' does not exist.", file=sys.stderr)
        raise typer.Exit(code=1)

    metadata: dict[str, Any] = read_metadata(metadata_toml)
    blast_defaults: dict[str, float ] = metadata["defaults"]["blast"]
    libraries: dict[str, Any] = metadata["libraries"]
    evalue: float = blast_defaults["evalue"]
    coverage: float = blast_defaults["coverage"]

    # Submit blasts in parallel
    vista_result: dict[str, Any] = dict()
    library_names = list(libraries.keys())
    workers, blast_threads = search_parallelism(cpus, len(library_names))
    search_func = partial(
        library_search,
        libraries=libraries,
        data_path=resources_path,
        evalue=evalue,
        coverage=coverage,
        query_fasta=query_fasta,
        num_threads=blast_threads,
    )
    with Pool(processes=workers) as pool:
        for library, result in pool.map(search_func, library_names):
            vista_result = vista_result | result

    print(json.dumps(vista_result, default=lambda x: x.__dict__), file=sys.stdout)


@app.command()
def build(
    data_path: Annotated[
        Path | None,
        typer.Option(
            "-d",
            "--data-path",
            help="The location of the input sequences and metadata. Defaults to './src/vista/config' or package resources.",
            file_okay=False,
            dir_okay=True,
            readable=True,
        ),
    ] = None,
    resource_dir: Annotated[
        Path | None,
        typer.Option(
            "-o",
            "--out-dir",
            help="The location of the BLAST databases. Defaults to './src/vista/resources' or package resources.",
            file_okay=False,
            dir_okay=True,
            readable=True,
        ),
    ] = None,
) -> None:
    """Generates the required BLAST databases."""
    if data_path is None:
        data_path = get_data_path()

    if resource_dir is None:
        resource_dir = get_resources_dir()

    if not data_path.is_dir():
        print(f"Error: data_path '{data_path}' does not exist.", file=sys.stderr)
        raise typer.Exit(code=1)

    if not resource_dir.exists():
        resource_dir.mkdir(parents=True, exist_ok=True)
    elif not resource_dir.is_dir():
        print(f"Error: output path '{resource_dir}' is not a directory.", file=sys.stderr)
        raise typer.Exit(code=1)

    metadata: dict[str, Any] = read_metadata(data_path / "metadata.toml")
    libraries = metadata["libraries"]

    # First, build simple marker databases
    for name in libraries.keys():
        if name != VIRULENCE_SETS:
            build_blastdb(
                resource_dir,
                name,
            )

    # Then build a DB for the virulence gene clusters
    virulence_genes: set[str] = set()
    for cluster in libraries[VIRULENCE_SETS]:
        virulence_genes.update(set(libraries[VIRULENCE_SETS][cluster]["genes"]))
    build_blastdb(resource_dir, VIRULENCE_SETS)


def get_data_path() -> Path:
    # If running from source, use local paths, otherwise use package resources
    if Path("pyproject.toml").exists():
        data_path = Path("src/vista/config")
    else:
        data_path = package_resource_path("config")
    return data_path


def get_resources_dir() -> Path:
    if Path("pyproject.toml").exists():
        resource_dir = Path("src/vista/resources")
    else:
        resource_dir = package_resource_path("resources")
    return resource_dir


if __name__ == "__main__":
    app()
