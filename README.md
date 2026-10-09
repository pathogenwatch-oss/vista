# Vista

[Change log](CHANGELOG.md) · [Evidence](EVIDENCE.md)

## Table of Contents

- [About](#about)
- [How to use](#how-to-use)
- [Installation](#installation)
- [Usage](#usage)
- [Example output](#example-output)
- [Acknowledgements](#acknowledgements)
- [Contributors](#contributors)
- [Licensing](#licensing)

## About

Vista is a database and genome assembly FASTA search tool for identifying _Vibrio cholerae_ serotypes,
along with identifying virulence genes and clusters, including ctx toxin genes.

This tool is currently under development by the [CGPS](https://www.pathogensurveillance.net/). Please open an issue or
contact us via [email](mailto:pathogenwatch@cgps.group) if you would link to know more or contribute.

## How to use

Vista takes a DNA sequence FASTA file as input and outputs a JSON format result to STDOUT.
It is recommended to install either install it as a python package, as a Docker image, or to run it directly with `uv`or
`pixi`. Vista has both `search` and `build` commands - most users will only need the search command.

```terminaloutput
 Usage: vista search [OPTIONS] QUERY_FASTA                                                                       
                                                                                                                 
╭─ Arguments ───────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *    query_fasta      FILE  The path to the query FASTA file. [required]                                      │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --metadata-toml  -m      FILE       The path to the metadata TOML file. Defaults to                           │
│                                     './src/vista/config/metadata.toml' or package resources.                  │
│ --data-path      -d      DIRECTORY  The location of the BLAST databases. Defaults to './src/vista/resources'  │
│                                     or package resources.                                                     │
│ --cpus           -c      INTEGER    The number of processes to run simultaneously. Defaults to the number of  │
│                                     CPUs.                                                                     │
│                                     [default: 8]                                                              │
│ --help                              Show this message and exit.                                               │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────
```

## Installation

First, clone this git repository. Note that a pre-compiled database is provided for easy installation and should work
for most users. If you need to rebuild the databases, see the [Building the databases](#building-the-database) section.

- `pixi` ensures a compatible version of BLAST will be installed along with python and all other required packages.
- `uv` will automatically install python and packages.
- `pip` Provided you have compatible python and BLAST versions installed to your system, vista can be installed as a
  system executable using pip.
- `Docker` can be used to create a portable versioned container.

```bash
git clone --depth 1 {repository}
cd vista
```

Then follow the most appropriate instructions for installation or running.

### Python/Pip

Running vista directly requires also installing required python packages, so it is recommended to install it using pip.

```bash
pip install . --no-cache-dir
cd /to/another/dir
vista search --help
```

### Pixi

Vista can be run with zero installation using pixi:

```bash
cd vista
pixi run vista search --help
```

To use it within a conda environment:

```bash
pixi install
pixi shell
vista search --help
```

### uv

You will need to install blastn separately. `uv` can be used to install `vista` as a module as well.

```bash
cd vista
uv run vista search --help
# From another directory
uv run --project /path/to/vista/repo/ vista search --help
```

### Docker

```bash
cd vista
docker build --rm -t vista .
cd ~/my_fasta_dir
docker run --rm -v $PWD:/fastas vista /fastas/my_vibrio_genome.fasta > result.json
```

### requirements.txt

A [requirements.txt](/requirements.txt) file is also provided to support alternate methods of installing the vista
module.

## Usage

### Building the database

See options with `vista build --help`.

```terminaloutput
❯ uv run vista build --help
                                                                                                                                                                                                                                   
 Usage: vista build [OPTIONS]                                                                                                                                                                                                      
                                                                                                                                                                                                                                   
 Generates the required BLAST databases.                                                                                                                                                                                           
                                                                                                                                                                                                                                   
╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --data-path  -d      DIRECTORY  The location of the input sequences and metadata. Defaults to './src/vista/config' or package resources.                                                                                        │
│ --out-dir    -o      DIRECTORY  The location of the BLAST databases. Defaults to './src/vista/resources' or package resources.                                                                                                  │
│ --help                          Show this message and exit.                                                                                                                                                                     │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

## Output description

### Output field descriptions

- `serogroup`: The predicted _Vibrio cholerae_ serogroup (for example, `O1`, `O139`, or
  `non-O1/O139`).
- `serogroupMarkers`: The list of searched markers, showing their names, the associated serogroup and the list of matches to the
  query genome.
- `virulenceGenes`: The list of virulence genes, with the following fields:
    - name: The name of the virulence gene.
    - type: The functional category of the gene (e.g. "Toxin", "Adhesion").
    - verified: Whether virulence impact has been experimentally verified in a mammalian host.
    - status: The presence status (`Present`, `Incomplete`, or `Not found`).
    - matches: A list of alignments found for that gene. Each match includes location details, identity, and whether the
      match is complete or disrupted.
- `virulenceClusters`: A list of objects representing virulence-associated gene clusters.
    - name: The name of the cluster.
    - genes: A list of all genes that are members of this cluster.
    - verified: Whether virulence impact has been experimentally verified in a mammalian host for at least one member of
      the cluster.
    - present/missing/incomplete: Lists of genes in the cluster categorised by their presence status.
    - status: An overall status for the cluster (`Present`, `Incomplete`, or `Not found`) based on the presence of its
      member genes.

### Example output

```json
{
  "virulenceGenes": [
    {
      "name": "ctxA",
      "type": "Toxin",
      "verified": "Yes",
      "status": "Present",
      "matches": [
        {
          "queryId": "ctxA",
          "contigId": "contig_1",
          "queryStart": 3385,
          "queryEnd": 4161,
          "refStart": 1,
          "refEnd": 777,
          "frame": 1,
          "isForward": true,
          "isComplete": true,
          "isDisrupted": false,
          "isExact": true,
          "identity": 100.0
        }
      ]
    }
  ],
  "virulenceClusters": [
    {
      "name": "Cqs quorum sensing cluster",
      "type": "Quorum sensing",
      "verified": "No",
      "genes": [
        "cqsS",
        "cqsA"
      ],
      "id": "cqs",
      "matches": {
        "cqsS": {
          "status": "Present",
          "matches": [
            {
              "queryId": "cqsS",
              "contigId": "contig_2",
              "queryStart": 1,
              "queryEnd": 2000,
              "refStart": 1,
              "refEnd": 2000,
              "frame": 1,
              "isForward": true,
              "isComplete": true,
              "isDisrupted": false,
              "isExact": true,
              "identity": 100.0
            }
          ]
        },
        "cqsA": {
          "status": "Not found",
          "matches": []
        }
      },
      "present": [
        "cqsS"
      ],
      "missing": [
        "cqsA"
      ],
      "incomplete": [],
      "status": "Incomplete"
    }
  ],
  "serogroup": "O1",
  "serogroupMarkers": [
    {
      "name": "rfbV",
      "type": "O1",
      "matches": [
        {
          "queryId": "rfbV",
          "contigId": "contig_3",
          "queryStart": 3385,
          "queryEnd": 4161,
          "refStart": 1,
          "refEnd": 777,
          "frame": 1,
          "isForward": true,
          "isComplete": true,
          "isDisrupted": false,
          "isExact": true,
          "identity": 100.0
        }
      ]
    },
    {
      "name": "wbfZ",
      "type": "O139",
      "matches": []
    }
  ]
}
```

## Acknowledgements

Originally developed by Corin Yeats and Sina Beier as part of the Vibriowatch project between
the [Centre for Pathogen Genome Surveillance](https://pathogensurveillance.net/), [Big Data Institute, Oxford](https://www.bdi.ox.ac.uk/)
and Nick Thompson's team at the [Wellcome Sanger Institute](https://www.sanger.ac.uk/). We would like to acknowledge the
support of our hosting institutes.

## Contributors

- Corin Yeats
- Sina Beier
- Avril Coghlan
- Nick Thompson
- David Aanensen

## Licensing

See [LICENCE](LICENSE).
