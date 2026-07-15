# syntax=docker/dockerfile:1.7

# Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
# Distributed under the terms of the Apache License, Version 2.0.

FROM rust:1-slim-bookworm AS rustbuild
RUN apt-get update && apt-get install -y --no-install-recommends build-essential pkg-config \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY modules/external/longread/Cargo.toml modules/external/longread/Cargo.lock ./
COPY modules/external/longread/src ./src
RUN cargo build --release --locked && strip target/release/longread
RUN cargo install xloci --version 0.0.4 --root /opt/xloci --locked

FROM mambaorg/micromamba:1.5.8-bookworm-slim
LABEL org.opencontainers.image.title="longread-all" \
      org.opencontainers.image.description="All longread pipeline tools in one image" \
      org.opencontainers.image.licenses="MIT"
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        pbsim3=3.0.4 pbccs=6.4.0 isoseq=4.0.0 lima samtools=1.21 procps-ng \
    && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:$PATH
COPY --from=rustbuild /src/target/release/longread /opt/conda/bin/longread
COPY --from=rustbuild /opt/xloci/bin/xloci /opt/conda/bin/xloci
COPY assets/pbsim_models /opt/pbsim_models
CMD ["bash"]
