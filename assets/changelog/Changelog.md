# Changelog

## v0.0.1 — 2026-07-15

### Summary

Initial public release of the longread pipeline — a reproducible, containerized Nextflow DSL2 workflow that simulates PacBio Iso-Seq reads from a BED12 transcript annotation and a reference genome.

### Added

#### Rust transcript and expression engine (`modules/longread/`)

- **CLI subcommands**: `prepare`, `validate`, `pbsim-input`, `split`, `rg`, `check`.
- **BED12 validation**: enforces unique transcript names, strand, sorted non-overlapping blocks, chromosome bounds, consistent gene-to-chromosome mapping, and no duplicate exon structures.
- **Synthetic isoform generation**: five local event types — alternative donor, alternative acceptor, exon skipping, intron retention, and five-prime truncation.
- **Fusion construction**: generates two-gene, same-chromosome, same-strand fusions with a configurable intronic gap. Outputs valid BED12 records and preserves 5'-to-3' partner order after strand-aware extraction.
- **Structural deduplication**: a global structural key (chromosome, strand, ordered exon starts and ends) prevents duplicate transcripts across originals, generated splice variants, truncations, and fusions.
- **Gene expression model**: GMM-derived raw weights scaled and normalized to a target molecule count using deterministic largest-remainder allocation.
- **Isoform allocation**: Zipf-like weights (`r^(-alpha)`) are assigned per-gene, normalized, and rounded to integer counts that sum exactly to the parent gene's total.
- **Deterministic parallelism**: per-gene and per-fusion seeds derived from `hash(global_seed, id)` via stable hashing. Byte-identical output across any thread count.
- **Output files**: `isoforms.bed`, `transcript_gene.tsv`, `gene_depth.tsv`, `isoform_depth.tsv`, `manifest.tsv`, `stats.json`.
- **PBSIM3 transcript input builder**: constructs the four-column transcript table from extracted sequences and isoform depths, omitting zero-count isoforms.
- **Chunking**: greedy largest-first bin packing for balanced PBSIM3 workload distribution.
- **BAM normalization**: rewrites PBSIM3 per-chunk BAMs into one synthetic PacBio movie with globally unique ZMWs and a specification-compliant SUBREAD read group.

#### Nextflow DSL2 pipeline (`src/`)

- **Entry point**: `main.nf` with `-params-file` and profile support.
- **Workflows**: `longread.nf` orchestrates the full data path.
- **Subworkflows**:
  - `prepare_transcriptome.nf` — validate, generate, extract sequences.
  - `simulate_subreads.nf` — chunk, run PBSIM3, merge subreads.
  - `process_isoseq.nf` — optional CCS and Iso-Seq clustering.
- **Local modules** (18 total): `chromsize`, `longread/prepare`, `longread/pbsim3`, `longread/split`, `xloci/exon`, `pbsim3`, `pbccs`, `merge_subreads`, `merge_ccs`, `validate_bam`, `pacbio/normalize_rg`, `pacbio/validate_ccs_chunks`, `bamtofa`, `isoseq/cluster2`, `publish`, `pbtk/pbindex`, `pbtk/pbmerge`.
- **Container strategy**: module-specific Dockerfiles and a whole-pipeline OCI image.
- **Profiles**: `local`, `docker`, `apptainer/singularity`, `slurm`, `test`.
- **Configuration**: `nextflow.config` with process-scoped resources, `nextflow_schema.json` for parameter validation.

#### Assets

- **Dockerfile**: whole-pipeline container image definition (`assets/docker/Dockerfile`).
- **PBSIM3 error models**: pre-trained `QSHMM-RSII`, `ERRHMM-SEQUEL`, and `ERRHMM-ONT` models.
- **Pipeline diagram**: Mermaid graph of the module graph (`assets/pipeline/longread.mermaid`).
- **SLURM runner script**: `assets/scripts/longread.sh` for cluster array-job submission.
- **Hiller Lab logo**: project branding (`assets/figures/hillerlab.png`).

#### Documentation

- `README.md` — project overview, usage instructions, output directory layout, configuration reference.
- `docs/usage.md` — detailed CLI and pipeline parameter documentation.
- `docs/decisions.md` — implementation notes and architectural decisions.
- `modules/longread/README.md` — crate-specific documentation.

#### Testing

- **Rust unit tests**: validation, prepare end-to-end, and PacBio-level (`modules/longread/tests/`).
- **Property tests**: random valid exon structures tested against block ordering, event correctness, deterministic output, and expression conservation.
- **Test data**: `mini.bed`, `mini.transcript_gene.tsv`, `genome.fa` for smoke tests and CI.

#### Infrastructure

- `LICENSE` — MIT license.
- `bin/.gitkeep` — placeholder for the compiled `longread` binary.
- `.gitignore` — Rust build artifacts, Nextflow work directories, and temporary files.

### Technical notes

- The Rust crate is pinned at version `0.0.6` (`edition = "2021"`, minimum Rust `1.81`).
- Key dependencies: `genepred` for BED parsing, `clap` for CLI, `rayon` for parallelism, `rand`/`rand_chacha`/`rand_distr` for deterministic RNG, `serde_json` for structured output, `noodles` for BAM I/O.
- All external tools are pinned to specific versions and invoked through container images built from pinned base-image digests.
