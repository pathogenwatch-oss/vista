# Vista

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
along with identifying virulence genes and clusters. It also provides ctxB allele assignments and copy numbers.

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

First clone this git repository. Note that a pre-compiled database is provided for easy installation and should work for
most users. If you need to rebuild the databases see the [Building the databases](#building-the-database) section.

- `pixi` ensures a compatible version of BLAST will be installed along with python and all other required packages.
- `uv` will automatically install python and packages.
- `pip` Provided you have compatible python and BLAST versions installed to your system vista can be installed as system
  executable using pip.
- `Docker` can be used to create a portable versioned container.

```bash
git clone --depth 1 {repository}
cd vista
```

Then follow the most appropriate instructions for installation or running.

### Python/Pip

Running vista directly requires also installing required python packages so it is recommended to install it using pip.

```bash
pip install . --no-cache-dir
cd /to/another/dir
vista search --help
```

### Pixi

Vista can be run with zero installation other than pixi with:

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
The output JSON contains the following main keys:
- `serogroup`: The predicted Vibrio cholerae serogroup (e.g., "O1", "O139", "Non-O1/O139").
- `serogroupMarkers`: A list of objects, one for each serogroup marker gene searched. Each object includes the gene name, the serogroup it's associated with, and a list of any matches found in the query genome.
- `virulenceGenes`: A list of objects, one for each virulence gene found.
  - name: The name of the virulence gene.
  - type: The functional category of the gene (e.g., "Toxin", "Adhesion").
  - status: The presence status ("Present", "Incomplete", or "Absent").
  - matches: A list of alignments found for that gene. Each match includes location details, identity, and whether the match is complete or disrupted.
- `virulenceClusters`: A list of objects representing virulence-associated gene clusters (e.g., TCP cluster).
  - name: The name of the cluster.
  - genes: A list of all genes that are members of this cluster.
  - present/missing/incomplete: Lists of genes in the cluster categorized by their presence status.
  - status: An overall status for the cluster ("Present", "Incomplete", or "Absent") based on the presence of its member genes.


### Example output

```json
{
  "virulenceGenes": [
    {
      "name": "ctxA",
      "type": "Toxin",
      "status": "Present",
      "matches": [
        {
          "contigId": ".CNRVC970056_CATTTT_L002.23",
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
      "name": "TCP cluster",
      "type": "colonisation",
      "genes": [
        "tcpA",
        "tcpB",
        "tcpC",
        "tcpD",
        "tcpE",
        "tcpF",
        "tcpH",
        "tcpI",
        "tcpJ",
        "tcpN",
        "tcpQ",
        "tcpR",
        "tcpS",
        "tcpT"
      ],
      "id": "tcp",
      "matches": {},
      "present": [
        "tcpA",
        "tcpB",
        "tcpC",
        "tcpD",
        "tcpE",
        "tcpF",
        "tcpH",
        "tcpI",
        "tcpJ",
        "tcpN",
        "tcpQ",
        "tcpR",
        "tcpS",
        "tcpT"
      ],
      "missing": [],
      "incomplete": [],
      "status": "Present"
    }
  ],
  "serogroup": "O1",
  "serogroupMarkers": [
    {
      "gene": "rfbV",
      "name": "O1",
      "matches": [
        "..."
      ]
    },
    {
      "gene": "wbfZ",
      "name": "O139",
      "matches": [
      ]
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
