# Auto-generated per-module image for nf-xen (clustering).
# Built from the module's conda environment — keep the two in sync.
FROM mambaorg/micromamba:1.5.10-noble

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH
USER root
WORKDIR /work
