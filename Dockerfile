FROM ubuntu:22.04 as blast

ARG BLAST_VERSION=2.9.0

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update \
      && apt install -y -q apt-transport-https software-properties-common curl \
      && rm -rf /var/lib/apt/lists/*

ENV BLAST_VERSION ${BLAST_VERSION}

RUN echo "Running ${BLAST_VERSION}"

RUN mkdir -p /tmp/blast \
    && mkdir /opt/blast \
    && curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/"${BLAST_VERSION}"/ncbi-blast-"${BLAST_VERSION}"+-x64-linux.tar.gz \
    | tar -zxC /tmp/blast --strip-components=1 \
    && cd /tmp/blast/bin \
  # The below line works for the latest versions of BLAST
  #  && cd /tmp/blast/ncbi-blast-"${BLAST_VERSION}"+/bin \
    && mv blastn makeblastdb /opt/blast/ \
    && cd .. \
    && rm -rf /tmp/blast

FROM python:3.10 as builder

COPY --from=blast /opt/blast/makeblastdb /opt/blast/makeblastdb

COPY data /data

ENV PATH /opt/blast:$PATH

RUN cd /data && \
      makeblastdb -in references.fasta -out references -dbtype nucl -parse_seqids

FROM ubuntu:22.04

RUN apt update \
    && apt install -y -q --no-install-recommends python3 python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install biopython toml \
    && pip3 cache purge \
    && apt -y -q remove python3-pip

COPY --from=builder /data /data

COPY --from=blast /opt/blast/blastn /opt/blast/blastn

COPY vista /vista

COPY vista.py /

ENV PATH /opt/blast:$PATH

ENTRYPOINT ["python3", "/vista.py"]