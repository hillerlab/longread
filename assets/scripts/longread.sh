#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# Example SLURM submission script for longread (nf-core pipeline)
#
# Each array task runs one reference × query genome pair. Nextflow itself is the
# "main job" — it submits all compute work as child SLURM jobs and only needs
# a small memory footprint.
#
# MANIFEST FILE FORMAT (species_list)
# ─────────────────────────────────────
# A tab-separated file with one pair per line, no header:
#
#   <query_prefix>  <query_bed>  <query_genome>  <query_gene_transcript_map>
#
# Example:
#   mm39   /data/genomes/mm39.bed   /data/genomes/mm39.2bit   /data/genomes/mm39.gene_map.tsv
#
# Genome files can be FASTA (.fa / .fasta) or 2bit (.2bit).
# Paths must be absolute. 
#
# USAGE
# ─────
# Edit the four path variables below, then submit with:
#   sbatch --array=1-<N> longread.sh
# where <N> is the number of lines in your manifest file.
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=makeChains
#SBATCH --array=1-10        # set upper bound to number of lines in species_list
#SBATCH -t 2-0
#SBATCH --output=/path/to/logs/%A.%a.out
#SBATCH --error=/path/to/logs/%A.%a.err
#SBATCH --mem=20G           # memory for the Nextflow process itself (not compute jobs)
#SBATCH -p public
#SBATCH -q public

# ── Load required modules (adjust to your cluster's module system) ────────────
module load nextflow
module load openjdk

# ── Environment ───────────────────────────────────────────────────────────────
export SLURM_SKIP_EPILOG=1

# Directory where Apptainer caches pulled container images
export NXF_APPTAINER_CACHEDIR=/scratch/$USER/longread/apptainer

# Give Nextflow's JVM enough heap for large runs (thousands of jobs)
export NXF_OPTS="-Xms4g -Xmx16g"

# ── Paths — edit these ────────────────────────────────────────────────────────
species_list="/path/to/manifest.tsv"   # tab-separated manifest (see format above)
working_dir="/path/to/output"          # one subdirectory per pair will be created here
pipeline_dir="/path/to/longread"  # cloned pipeline repo

# ── Parse manifest line for this array task ───────────────────────────────────
pair=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$species_list")
query_prefix=$(echo "$pair" | cut -f1)
query_bed=$(echo "$pair"    | cut -f2)
query_genome=$(echo "$pair" | cut -f3)
query_gene_transcript_map=$(echo "$pair" | cut -f4)

if [[ -z "$query_prefix" || -z "$query_bed" || -z "$query_genome" || -z "$query_gene_transcript_map" ]]; then
    echo "ERROR: could not parse line ${SLURM_ARRAY_TASK_ID} of ${species_list}" >&2
    exit 1
fi

# ── Per-pair working directory ─────────────────────────────────────────────────
pair_dir="${working_dir}/${query_prefix}_isoseq_reads"
mkdir -p "${pair_dir}/logs"

# ── Write params.json for this pair ───────────────────────────────────────────
# Scientific parameters go here; infrastructure stays in nextflow.config.
cat > "${pair_dir}/params.json" <<EOF
{
    "prefix":   "${query_prefix}",
    "bed":    "${query_bed}",
    "sequence": "${query_genome}",
    "transcript_gene": "${query_gene_transcript_map}",
    "outdir":        "${pair_dir}/results",

    "use_container": true,

    "seed": 42,

    "new_isoforms_per_gene": null,
    "mean_new_isoforms_per_gene": 1.0,
    "max_new_isoforms_per_gene": 5,
    "event_weights": "donor=1,acceptor=1,skip=1,retention=1,truncation=1",
    "max_event_attempts": 100,
    "minimum_transcript_length": 200,

    "fusion_count": 0,
    "min_fusion_intron": 1,
    "max_fusion_distance": null,
    "fusion_expression_scale": 1.0,

    "total_molecules": 100000,
    "alpha": 4.0,

    "pbsim_method": "errhmm",
    "pbsim_model": "${projectDir}/../assets/pbsim_models/ERRHMM-SEQUEL.model",
    "pass_count": 10,
    "pbsim_chunks": 1000,
    "pbsim_machine": "SEQUEL",

    "do_ccs": true,
    "do_isoseq_cluster": true,
    "ccs_chunk": 1000,
    "isoseq_min_rq": 0.95
}
EOF

# cd into pair_dir so each run's .nextflow.log is saved there
cd "$pair_dir"

nextflow run "${pipeline_dir}/main.nf" \
    -params-file "${pair_dir}/params.json" \
    -profile     apptainer,slurm \
    -w           "${pair_dir}/work"
