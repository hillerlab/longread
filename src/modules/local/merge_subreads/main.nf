/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
 * Merge per-chunk subread BAMs into one ${prefix}.subreads.bam (CLAUDE.md §5.5).
 */

process MERGE_SUBREADS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23--h96c455f_0' :
        'biocontainers/samtools:1.23--h96c455f_0' }"

    input:
    tuple val(meta), path(bams)
    tuple val(meta1), path(mafs)

    output:
    tuple val(meta), path("${meta.id}.subreads.bam"), emit: bam
    tuple val(meta), path("${meta.id}.subreads.maf"), emit: maf
    path "versions.yml",                              emit: versions

    script:
    """
    samtools merge -@ ${task.cpus} -f ${meta.id}.subreads.bam ${bams}
    samtools quickcheck -u ${meta.id}.subreads.bam

    cat ${mafs} > ${meta.id}.subreads.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.subreads.bam
    touch ${prefix}.subreads.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
