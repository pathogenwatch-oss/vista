FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update \
    && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY vista /vista

RUN mkdir /data

COPY data/references.fasta /data/

COPY data/metadata.toml /data/

RUN cd /data && makeblastdb -in references.fasta -out references -dbtype nucl -parse_seqids

RUN pip3 install biopython toml

COPY vista.py /

#CMD ["/bin/bash"]
ENTRYPOINT ["python3", "/vista.py"]