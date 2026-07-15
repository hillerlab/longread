# longread

A reproducible, containerized **Nextflow DSL2** pipeline that simulates PacBio Iso-Seq data from a
BED12 transcript annotation, a transcript→gene mapping, and a reference genome.

```text
BED12  +  transcript_gene.tsv  +  genome(.fa/.fa.gz/.2bit)
                       │
                       ▼
   gene-aware alternative splicing + 5′ truncation + two-gene fusions   (longread prepare)
                       │
                       ▼
        exact gene-level GMM + within-gene isoform allocation
                       │
                       ▼
        full spliced-transcript extraction                              (xloci)
                       │
                       ▼
        PBSIM3 --strategy trans  →  merged subreads BAM                 (pbsim3 + samtools)
                       │
              optional CCS (HiFi)                                       (ccs)
                       │
              optional Iso-Seq clustering                               (isoseq cluster2)
```

The engine is a native Rust tool (`modules/external/longread/`, `bin/longread`) built on the
`genepred` BED parser; sequence extraction uses `xloci`. See `CLAUDE.md` for the full design and
`docs/decisions.md` for implementation notes.

## Quick start

```bash
# From the repository root:
nextflow run src/main.nf -profile test,local        # tools from PATH
nextflow run src/main.nf -profile test,docker        # per-module containers
nextflow run src/main.nf -profile test,docker --do_ccs --do_isoseq_cluster
```

Real run:

```bash
nextflow run src/main.nf -profile docker \
    --bed transcripts.bed \
    --transcript_gene transcript_gene.tsv \
    --sequence genome.2bit \
    --chrom_sizes genome.chrom.sizes \
    --prefix mysim \
    --total_molecules 1000000 \
    --fusion_count 100 \
    --pbsim_chunks 16 \
    --pass_count 10 \
    --do_ccs --do_isoseq_cluster \
    --outdir results
```

## Inputs

| Parameter | Description |
|-----------|-------------|
| `--bed` | BED12 transcript annotation (required) |
| `--transcript_gene` | `TRANSCRIPT_ID⇥GENE_ID` mapping (required) |
| `--sequence` | reference genome `.fa` / `.fa.gz` / `.2bit` (required) |
| `--chrom_sizes` | chromosome sizes; enables coordinate-bounds validation (optional) |

Key knobs: `--new_isoforms_per_gene` **or** `--mean_new_isoforms_per_gene`/`--max_new_isoforms_per_gene`,
`--event_weights`, `--fusion_count`, `--total_molecules`, `--alpha`, `--pbsim_method`, `--pass_count`,
`--pbsim_chunks`, `--do_ccs`, `--do_isoseq_cluster`, `--seed`. Full list: `nextflow run src/main.nf --help`.

## Outputs (`--outdir`)

| File | Description |
|------|-------------|
| `${prefix}.subreads.bam` | merged PBSIM3 subreads (always) |
| `${prefix}.ccs.bam` | merged CCS/HiFi reads (`--do_ccs`) |
| `${prefix}.clustered.bam` | Iso-Seq transcript clusters (`--do_isoseq_cluster`) |
| `prepare/${prefix}.{isoforms.bed,isoform_depth.tsv,gene_depth.tsv,manifest.tsv,stats.json,transcript_gene.tsv}` | transcript inventory & expression |
| `chunks/chunks.tsv` | chunk manifest |
| `validation/*.validation.txt`, `pipeline_info/*` | BAM validation, reports, `software_versions.yml` |

## Profiles

`local` (tools from PATH), `docker`, `singularity`/`apptainer`, `slurm` (combine, e.g. `-profile slurm,singularity`),
`test` (bundled fixtures under `tests/data/`). Resources are clamped by `--max_cpus/--max_memory/--max_time`.

## Testing

```bash
# Rust engine (Parts 1-2): unit + property + integration tests
cd modules/external/longread && cargo test && cargo clippy --all-targets -- -D warnings && cargo fmt --check

# Tool integration (Part 2, §5.8): xloci → pbsim → split → PBSIM3 → merge → CCS → isoseq
tests/integration/part2_pipeline.sh

# Pipeline (Part 3, §7 Test F)
nextflow run src/main.nf -profile test,local
```

## Containers

Each external tool has a dedicated image (`modules/external/*/Dockerfile`); an all-in-one image is at
`assets/docker/Dockerfile`. The `docker`/`singularity` profiles default to pinned biocontainers for the
bioconda tools and `hillerlab/{longread,xloci}` images for the Rust tools. Pin to digests per release
(`CLAUDE.md` §6.4).

## License

MIT (see `LICENSE`). Third-party tool licenses: `THIRD_PARTY_LICENSES.md`.
