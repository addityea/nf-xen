# Per-module image for nf-xen (celltyping_singler).
FROM mambaorg/micromamba:1.5.10-noble

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH
USER root
WORKDIR /work
