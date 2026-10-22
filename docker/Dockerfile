# nf-xen Python stack: scanpy + squidpy + celltypist + report helpers.
# Build:  docker build -t saditya88/nf-xen:python-1.0 docker/
FROM python:3.11-slim

ENV PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libhdf5-dev libgeos-dev libproj-dev pkg-config build-essential procps \
    && rm -rf /var/lib/apt/lists/*

# Pinned major stack for reproducibility
RUN pip install \
        anndata==0.11.4 \
        pandas==2.3.3 \
        scanpy==1.11.2 \
        squidpy==1.6.5 \
        celltypist==1.7.1 \
        leidenalg==0.10.2 \
        igraph==0.11.8 \
        scikit-learn \
        pyclustree \
        matplotlib \
        seaborn \
        pypdf \
        pyarrow

WORKDIR /work
