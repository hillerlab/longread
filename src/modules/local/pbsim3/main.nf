/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
 * PBSIM3 transcript-mode read simulation (CLAUDE.md §5.3-§5.5).
 * One task per chunk; a unique --id-prefix keeps read (movie/ZMW) names distinct at merge time.
 */

process PBSIM3 {
    tag "${meta.id}:${chunk}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/hillerlab/pbsim3:latest' }"

    input:
    tuple val(meta), path(chunk)
    path model

    output:
    tuple val(meta), path("*.bam"),           emit: bam
    tuple val(meta), path("*.subreads.maf"),  emit: maf
    tuple val(meta), path("*.ref_map.tsv"),   emit: map
    path "versions.yml",                      emit: versions

    script:
    def idpfx  = "${meta.id}"
    def out_bam = "${meta.id}.subreads.bam"
    def out_sam = "${meta.id}.subreads.sam"
    def out_maf = "${meta.id}.subreads.maf"
    def ref_map = "${meta.id}.ref_map.tsv"

    // Derive a deterministic per-chunk seed so error streams differ between chunks.
    def seed = params.seed + (meta.id.toString().hashCode() & 0x7fffffff) % 1000000

    def mode = params.pbsim_mode == "trans" ? "trans" : "wgs"
    def input = params.pbsim_mode == "trans"
        ? "--transcript ${chunk}"
        : "--genome ${chunk}"
    """
    pbsim3 \\
        --strategy ${mode} \\
        --method ${params.pbsim_method} \\
        ${input} \\
        --${params.pbsim_method} ${model} \\
        --pass-num ${params.pass_count} \\
        --prefix ${idpfx} \\
        --id-prefix movie.${idpfx} \\
        --seed ${seed}
    
    # if mode == wgs, merge all local files into one local chunk:
    if [ ${mode} == "wgs" ]; then
      for sam in *.sam; do
          samtools view \\
              -b \\
              -o "\${sam%.sam}.bam" \\
              "\${sam}"

          rm "\${sam}"
      done

      # wgs mode writes one .maf/.ref per transcript: never expand them on argv.
      find . -maxdepth 1 -name '*.maf' -print0 | sort -z | xargs -0 -r cat > maf.tmp
      find . -maxdepth 1 -name '*.maf' -delete
      mv maf.tmp ${out_maf}

      # Map each split reference file to its FASTA entry.
      printf 'ref_file\\tfasta_entry\\n' > ${ref_map}

      for ref in *.ref; do
          awk -v file="\$(basename "\${ref}")" '
              /^>/ {
                  print file "\\t" substr(\$0, 2)
              }
          ' "\${ref}" >> ${ref_map}
      done

      find . -maxdepth 1 -name '*.ref' -delete
    else
      # PBSIM3 emits multipass reads as BAM or SAM depending on the build; normalise to BAM.
      if [ -f ${idpfx}.bam ]; then
          mv ${idpfx}.bam ${out_bam}
      elif [ -f ${idpfx}.sam ]; then
          samtools view -b -o ${out_bam} ${idpfx}.sam
      else
          echo "ERROR: PBSIM3 produced neither ${idpfx}.bam nor ${idpfx}.sam" >&2
          exit 1
      fi

      # trans mode writes <prefix>.maf whose reference names are the transcript IDs.
      mv ${idpfx}.maf ${out_maf}
      printf 'ref_file\\tfasta_entry\\n' > ${ref_map}
      awk -F'\\t' '{ print \$1 "\\t" \$1 }' ${chunk} >> ${ref_map}
    fi

    # PBSIM3 has no --version flag; the version is fixed by the pinned container image.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: 3.0.4
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch *.bam
    touch ${prefix}.subreads.maf
    touch ${prefix}.ref_map.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbsim3: 1.0.4
    END_VERSIONS
    """
}
