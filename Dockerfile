FROM ubuntu:22.04 as blast

ARG BLAST_VERSION=2.12.0

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update \
      && apt install -y -q apt-transport-https software-properties-common curl \
      && rm -rf /var/lib/apt/lists/*

ENV BLAST_VERSION ${BLAST_VERSION}

RUN echo "Using BLAST version: ${BLAST_VERSION}"

# The below line works for the latest versions of BLAST
#  && cd /tmp/blast/ncbi-blast-"${BLAST_VERSION}"+/bin \
RUN mkdir -p /tmp/blast \
    && mkdir /opt/blast \
    && curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/"${BLAST_VERSION}"/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz | tar -zxC /tmp/blast --strip-components=1 \
    && cd /tmp/blast/bin \
    && mv blastn makeblastdb /opt/blast/ \
    && cd .. \
    && rm -rf /tmp/blast

FROM ubuntu:22.04

COPY requirements.txt /

RUN apt update \
    && apt install -y -q --no-install-recommends python3 python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install -r /requirements.txt \
    && pip3 cache purge \
    && apt -y -q remove python3-pip

COPY --from=blast /opt/blast/makeblastdb /opt/blast/makeblastdb

COPY --from=blast /opt/blast/blastn /opt/blast/blastn

ENV PATH /opt/blast:$PATH

COPY data /data

COPY vista /vista

COPY vista.py /

RUN python3 /vista.py build

ENTRYPOINT ["python3", "/vista.py", "search"]