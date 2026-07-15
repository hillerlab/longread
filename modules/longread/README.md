# `longread` — transcript inventory & expression engine

The Rust engine behind Part 1 of the `longread` PacBio Iso-Seq simulation pipeline. It reads a
BED12 annotation plus a transcript→gene mapping, generates synthetic splice/truncation isoforms
and two-gene fusions, and assigns **exact** integer molecule counts at the gene and isoform level.

All randomness is deterministic in `(seed, gene_id)`, so a given input and seed produce
byte-identical output regardless of the thread count.

## Build

```bash
cargo build --release            # produces target/release/longread
cargo test                       # unit + property + integration tests
cargo clippy --all-targets --all-features -- -D warnings
cargo fmt --check
```

The pipeline exposes the compiled binary at `bin/longread`.

## Commands

### `longread prepare`

Builds the final transcript inventory and expression tables.

```bash
longread prepare \
    --bed transcripts.bed \
    --transcript-gene transcript_gene.tsv \
    --chrom-sizes genome.chrom.sizes \
    --output-prefix simulation \
    --mean-new-isoforms-per-gene 1.0 \
    --max-new-isoforms-per-gene 5 \
    --fusion-count 100 \
    --total-molecules 100000 \
    --alpha 4.0 \
    --seed 42 \
    --threads 16
```

Use `--new-isoforms-per-gene N` for a fixed count instead of the Poisson `mean`/`max` pair
(the two modes are mutually exclusive). `--event-weights donor=1,acceptor=1,skip=1,retention=1,truncation=1`
tunes event selection. `--require-exact-event-count` turns a shortfall into an error.

Outputs (spec §4.11):

| File | Contents |
|------|----------|
| `${prefix}.isoforms.bed` | BED12 of all transcripts (originals + generated + fusions) |
| `${prefix}.transcript_gene.tsv` | transcript→gene mapping (headerless) |
| `${prefix}.gene_depth.tsv` | `GENE_ID  GENE_DEPTH  RAW_WEIGHT  IS_FUSION_GENE` |
| `${prefix}.isoform_depth.tsv` | `TRANSCRIPT_ID  GENE_ID  ISOFORM_DEPTH  ANTISENSE_DEPTH  RELATIVE_WEIGHT` |
| `${prefix}.manifest.tsv` | per-transcript provenance |
| `${prefix}.stats.json` | run statistics |

Invariants guaranteed:

```text
sum(gene_depth)               == total_molecules
sum(isoform_depth for gene G) == gene_depth[G]
```

### `longread validate`

Validates the BED12 + mapping (+ optional chrom.sizes) and reports every problem at once.

```bash
longread validate --bed transcripts.bed --transcript-gene transcript_gene.tsv --chrom-sizes genome.chrom.sizes
```

### `longread pbsim`

Builds the PBSIM3 four-column transcript-mode input (`TRANSCRIPT_ID SENSE ANTISENSE SEQUENCE`) from
`xloci --as-tsv` sequences and `isoform_depth.tsv`. Zero-count transcripts are omitted; with
`--manifest` it also checks every fusion sequence equals its `5′ + 3′` partner concatenation.

```bash
longread pbsim --sequences sim.tsv --isoform-depth sim.isoform_depth.tsv \
    --manifest sim.manifest.tsv --output sim.pbsim.transcript.tsv
```

### `longread split`

Greedy largest-first chunking of the PBSIM3 transcript file into balanced chunks
(`work = length × sense × pass_count`), writing `chunk_NNNN.transcript.tsv` and `chunks.tsv`.

```bash
longread split --transcript sim.pbsim.transcript.tsv --pbsim-chunks 8 --pass-count 10 \
    --outdir chunks --prefix simulation
```

Each chunk's `CHUNK_PREFIX` must be passed to PBSIM3 as **`--id-prefix`** (not just `--prefix`) so
read names stay globally unique after merging.

### `longread rg`

Canonicalizes PBSIM3 subread BAMs as one synthetic PacBio movie: assigns globally unique ZMWs
(deterministic per-movie offsets), rewrites each record's QNAME and `zm` tag, and collapses all
records onto one specification-compliant `SUBREAD` read group (`ID = md5(lower(movie)//SUBREAD)[:8]`).
Emits one `*.normalized.bam` per input plus `zmw_map.tsv`. Runs in-process (noodles + libdeflate),
two parallel passes, deterministic across thread counts. ZMW allocation is keyed **per source BAM**
(each PBSIM3 file has its own restarted ZMW space), so movie-name collisions across files — as
produced by PBSIM3 `wgs` mode — are disambiguated rather than rejected.

`--bams` is a **list file** with one BAM path per line (avoids passing thousands of paths on argv):

```bash
ls -1 *.bam > bam.list
longread rg --movie movie.sample --bams bam.list --outdir . --zmw-map zmw_map.tsv --threads 16
```

`--max-open-files` (default 512) bounds simultaneously open file descriptors: work is parallelized
across BAMs with a pool of `min(threads, max_open_files/2)` workers, so the FD peak is independent
of the number of chunks.

### `longread check`

Fails before merge if CCS chunks overlap in PacBio's `movie/zmw` key space (a ZMW must occur in
exactly one chunk). Parallel scan + global de-duplication; writes a marker file on success.

```bash
ls -1 *.bam > bam.list
longread check --bams bam.list --out ccs_chunks.valid --threads 16
```

## Status

Parts 1 and 2 (this crate's subcommands) are complete and pass their acceptance gates
(`CLAUDE.md` §4.12, §5.8). The end-to-end tool integration (xloci → PBSIM3 → CCS → Iso-Seq) is
verified by `tests/integration/part2_pipeline.sh`. See `../../../docs/decisions.md` for notes.
