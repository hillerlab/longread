# longread engine — multi-stage build of the Rust binary.
# Build context: this directory (modules/external/longread/).
#   docker build -t hillerlab/longread:0.1.0 modules/external/longread
FROM rust:1-slim-bookworm AS build
WORKDIR /src
# gcc: the `libdeflate-sys` crate (fast bgzf codec used by the PacBio BAM subcommands)
# compiles bundled C sources with the `cc` crate.
RUN apt-get update && apt-get install -y --no-install-recommends pkg-config gcc \
    && rm -rf /var/lib/apt/lists/*
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release --locked && strip target/release/longread

FROM debian:bookworm-slim
LABEL org.opencontainers.image.title="longread" \
      org.opencontainers.image.description="BED12-native transcript & expression engine" \
      org.opencontainers.image.licenses="MIT"
# procps provides `ps`, required by Nextflow for task metrics.
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/longread /usr/local/bin/longread
USER 1000:1000
CMD ["longread", "--help"]
