# Vista

## WARNING

We have identified an issue with the selection of reference loci for biotyping in previous versions of Vista. These
results should not be used and are now excluded from the current version.

## About

A database and genome assembly search tool for identifying _Vibrio cholerae_ biotypes and serotypes, and identifying
virulence genes and clusters.

This tool is currently under development by the [CGPS](https://www.pathogensurveillance.net/). Please open an issue or
contact us via [email](mailto:pathogenwatch@cgps.group) if you would link to know more or contribute.

## How to use

Vista takes a DNA sequence FASTA file as input and outputs a JSON format result to STDOUT.
While it will run directly on the command line, we only support building and running Docker images.
To create a local or bespoke build, we recommend building a Docker image using the provided [Dockerfile](/Dockerfile).
However, it is straightforward to install and run locally using Python3 + BLAST.

### Running via Docker

```
git clone [vista repo]
cd vista
docker build --rm -t vista .
cd ~/my_fasta_dir
docker run --rm -v $PWD:/fastas vista /fastas/my_vibrio_genome.fasta > result.json
```

### Running directly

First install NCBI-BLAST (blastn & makeblastdb), python 3 and pip.

```
git clone [vista repo]
cd vista
pip3 install -r requirements.txt
python3 vista.py build
python3 vista.py search /path/to/my_vibrio_genome.fasta > result.json
```

## Example output

```
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
        },
        ... etc ...
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
            "matches": {...as above...},
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
        },
        ...etc...
    ],
    "serogroup": "O1",
    "serogroupMarkers": [
        {
            "gene": "rfbV",
            "name": "O1",
            "matches": [...as above...]
        },
        {
            "gene": "wbfZ",
            "name": "O139",
            "matches": [...as above...]
        }
    ],
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

See [LICENSE](LICENSE).
