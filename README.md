# Vista
## About
A database and genome assembly search tool for identifying _Vibrio cholerae_ biotypes and serotypes, and identifying virulence genes and clusters.

This tool is currently under development by the [CGPS](https://www.pathogensurveillance.net/). Please open an issue or contact us via [email](mailto:pathogenwatch@cgps.group) if you would link to know more or contribute.

## How to use
Vista takes a DNA sequence FASTA file as input and outputs a JSON format result to STDOUT. 
While it will directly on the command line, we only support running one of the official Vista Docker images. 
The current images can be found at [our GitLab registry](https://gitlab.com/cgps/vista/container_registry/893140). 
To create a local or bespoke build, we recommend building a Docker image using the provided [Dockerfile](/Dockerfile).
However, it is straightforward to install and run locally using Python3 + BLAST.

### Running via Docker
```
docker run --rm -v $PWD:/data/ registry.gitlab.com/cgps/vista:v0.0.6 my_vibrio_genome.fasta > result.json
```

### Running directly
```
python3 vista.py my_vibrio_genome.fasta $PWD > result.json
```

*NB* If the second argument is not provided, then it assumes the FASTA file is in /data.

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
    "biotype": "O1 El Tor",
    "biotypeMarkers": [
        {
            "name": "O1",
            "gene": "rfbV",
            "matches": [...as above...]
        },
        {
            "name": "O1 El Tor",
            "gene": "ctxB3",
            "matches": [...as above...]
        }
    ]
}
```

## Acknowledgements
Originally developed by Corin Yeats and Sina Beier as part of the Vibriowatch project between the Centre for Pathogen Genome Surveillance, Big Data Institute, Oxford and Nick Thompson's team at the Wellcome Sanger Institute. We would like to acknowledge the support of our hosting institutes.

## Contributors
 - Corin Yeats
 - Sina Beier
 - Avril Coghlan
 - Nick Thompson
 - David Aanensen

## Licensing
See [LICENSE](LICENSE).
