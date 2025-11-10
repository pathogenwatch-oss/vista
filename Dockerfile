FROM ghcr.io/astral-sh/uv:python3.13-bookworm-slim AS compile

ARG BLAST_VERSION=2.17.0

# Install build dependencies and BLAST in a single layer
RUN apt update && \
    apt install -y --no-install-recommends \
        curl \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /tmp/blast /opt/blast && \
    curl -fsSL "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz" | \
    tar -zxC /tmp/blast --strip-components=1 && \
    mv /tmp/blast/bin/blastn /tmp/blast/bin/makeblastdb /opt/blast/ && \
    rm -rf /tmp/blast

WORKDIR /vista

COPY src /vista/src
COPY pyproject.toml uv.lock LICENSE README.md /vista/

RUN uv build --wheel

FROM ghcr.io/astral-sh/uv:python3.13-bookworm-slim AS code

# Install runtime dependencies
RUN apt update && \
    apt install -y --no-install-recommends libgomp1 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /vista

COPY --from=compile /vista/LICENSE /vista/README.md /vista/
COPY --from=compile /vista/dist/vista*.whl /vista/dist/

RUN uv pip install --system --no-cache-dir /vista/dist/*.whl && \
    rm -rf /vista/dist

FROM code AS prod

COPY --from=compile /opt/blast/blastn /usr/local/bin/

ENTRYPOINT ["vista", "search"]