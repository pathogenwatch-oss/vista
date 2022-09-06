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

*NB* If the second option is not provided then it assumes the FASTA file is in /data.

## Example output


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
