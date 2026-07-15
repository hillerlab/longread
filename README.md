
<p align="center">
  <p align="center">
    <img width=200 align="center" src="./assets/figures/hillerlab.png" >
  </p>

  <span>
    <h1 align="center">
        longread
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/longread" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/longread?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> simulate PacBio Iso-Seq reads at scale </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://www.pacb.com/products-and-services/applications/rna-sequencing/">isoseq</a> .
        <a href="https://github.com/hillerlab/longread/blob/main/assets/pipeline/longread.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

<div align="center">

<pre style="font-size: 18px;">

REF: ══▓▓▓▓══░░░░░░══▓▓▓▓══░░░░░░══▓▓▓▓══

FL:  ══▓▓▓▓══░░░░░░══▓▓▓▓══▓▓▓▓══▓▓▓▓══

IR:  ══▓▓▓▓══▓▓▓▓▓▓▓▓▓▓▓▓══▓▓▓▓
</pre>

</div>

---

> [!IMPORTANT]
> - **Scaffold names**: Input genome is renamed. See []() for the specific formatting rules.
> - **Inputs accepted**: `.fasta`, `.2bit`, or `.gz`.
> - **Container image**: We offer a pre-built container image for the whole pipeline as well as individual modules. By default the pipeline runs with [ghcr.io/hillerlab/longread:latest](https://github.com/hillerlab/containers/pkgs/container/longread). Additional images can be found at [containers](https://github.com/hillerlab/containers) and nextflow modules at [core](https://github.com/hillerlab/core).

---

## Usage

> [!NOTE]
> Requirements: Nextflow ≥ 25.04.6, Docker or Apptainer, Java.

```bash
git clone https://github.com/hillerlab/longread.git
cd longread
```

Edit `params.json` (set `bed`, `sequence`, `transcript_gene`, `outdir`, `prefix`), then:

```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

Smoke test:
```bash
nextflow run main.nf -profile test,apptainer
```

> [!NOTE]
> You can also specify these options directly in `params.json`.

A helper sh script is provided to run the pipeline on a SLURM cluster. See details below.

<details>
<summary>Click to expand</summary>


Edit the path variables at the top of `assets/hpc/longread.sh` (cache dir, container image, manifest path), then submit:

```bash
sbatch --array=1-<N> do_longread.sh
```

Each array task spawns one Nextflow head job that submits all compute as child SLURM jobs.

REPEATMASKER run as SLURM job arrays. Partition routing, array sizes, and resource tiers are documented inline in `nextflow.config` — edit there to match your cluster.

</details>

---

## Output

```
results/
├── 01_CHROMSIZE/        chrom.sizes
├── 02_LONGREAD_PREPARE  *.tsv/*.bed/*.json
├── 03_XLOCI             *.fa
├── 04_PBSIM3
├────── BAM/             *.bam
├────── MAF/             *.maf
├────── MAP/             *.tsv
├── 05_MERGE_SUBREADS    *.bam
├── 06_PBCCS             *.bam
├── 07_ISOSEQ_CLUSTER    *.bam
├── 08_BAM_TO_FA         *.fa
└── PIPELINE_INFO/    timeline, trace, DAG, versions
```

---

## Where to edit

| File | What |
|------|------|
| `params.json` | Genome paths, alignment settings, checkpoints — per run |
| `nextflow.config` | Compute resources, profiles, container, SLURM — rarely |
