# Changelog

## v0.0.2 — 2026-10-09

### Summary

Scaling and robustness pass over the simulation path: PBSIM3 task sizing no longer doubles as a chunk count in wgs mode, chunks carry unique ids so prefixes and outputs never collide, and the modules that stage many inputs read them from the task directory instead of expanding them onto the command line.

### Added

#### Rust transcript and expression engine (`modules/longread/`, v0.0.7)

- **Molecule cap per transcript record**: `split` now breaks records heavier than 50,000 molecules (sense or antisense) into near-equal rows that keep the same transcript id, conserving the total counts exactly. This works around PBSIM3 3.0.4 trans mode segfaulting around ~60,000 molecules of a single transcript and additionally spreads heavy transcripts across chunks during bin packing.

#### Nextflow DSL2 pipeline (`src/`)

- **`pbsim_records_per_chunk` parameter** (default `1000`): number of FASTA records (transcripts) per PBSIM3 task in wgs mode. `pbsim_chunks` (default `10`) is now trans-mode only.
- **Unique chunk ids in trans mode**: each chunk gets `<sample>.chunk_N`, so PBSIM3 prefixes, derived seeds, chunk MAF names and subread BAM names never collide across chunks.
- **trans-mode `ref_map` output**: the PBSIM3 module now emits the chunk MAF and the `ref_map.tsv` (transcript id → fasta entry) that `merge_subreads` expects, matching the wgs code path.

### Changed

- Pipeline manifests bumped to `v0.0.2` (`nextflow.config`), engine crate to `0.0.7`.
- `test` profile: `pbsim_chunks = 2` for trans mode, `pbsim_records_per_chunk = 2` for wgs mode.
- Whole-pipeline image (`assets/docker/Dockerfile`) builds the engine from `modules/longread`, installs `pbtk 3.5.0` (required by the `pbtk/pbindex` and `pbtk/pbmerge` modules), symlinks bioconda's `pbsim` to `pbsim3` (the module invokes `pbsim3`, the package installs `pbsim`), and bumps `xloci` to `0.0.7`.

### Fixed

- **Staged inputs are listed from disk**: `isoseq/cluster2`, `merge_ccs`, `merge_subreads` and `pacbio/validate_ccs_chunks` build `bam.list`/`.fofn` with `find ... | sort` instead of interpolating staged BAMs into the process script, so the pipeline no longer grows command lines with chunk counts.
- **Per-transcript MAF handling in pbsim3**: wgs mode concatenates with `find -print0 | sort -z | xargs -0 cat` and deletes with `find -delete` instead of `cat *.maf` / `rm *.maf`.

### Technical notes

- `docs/usage.md` documents the `pbsim_mode`/`pbsim_chunks`/`pbsim_records_per_chunk` knob split.

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
