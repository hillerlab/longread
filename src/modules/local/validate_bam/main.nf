/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
 * Validate a BAM: samtools quickcheck (-u for unaligned PacBio BAMs) + header dump.
 * Fails the pipeline if the BAM is truncated or malformed (CLAUDE.md §5.5).
 */

process VALIDATE_BAM {
    tag "${meta.id}:${bam.name}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.23--h96c455f_0' :
        'biocontainers/samtools:1.23--h96c455f_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.baseName}.validation.txt"), emit: report
    path "versions.yml",                                     emit: versions

    script:
    """
    samtools quickcheck -u ${bam}
    {
        echo "file: ${bam}"
        echo "reads: \$(samtools view -c ${bam})"
        echo "quickcheck: OK"
        echo "--- header ---"
        samtools view -H ${bam}
    } > ${bam.baseName}.validation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.subreads.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
