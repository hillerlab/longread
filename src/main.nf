#!/usr/bin/env nextflow

/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    longread

    Simulate PacBio Iso-Seq reads at scale
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller
    GitHub:  https://github.com/hillerlab/longread
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS AND MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LONGREAD as SIMULATE } from './workflows/longread.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateFullRun() {
    def problems = []

    if (!params.bed)             { problems << 'missing required --bed' }
    if (!params.transcript_gene) { problems << 'missing required --transcript_gene' }
    if (!params.sequence)        { problems << 'missing required --sequence (genome)' }

    if (params.do_ccs && (params.pass_count as int) < 2) {
        problems << "do_ccs requires pass_count >= 2 (got ${params.pass_count})"
    }

    if (params.do_isoseq_cluster && !params.do_ccs) {
        problems << 'do_isoseq_cluster requires do_ccs = true'
    }

    if (!(params.pbsim_method in ['errhmm', 'qshmm'])) {
        problems << "pbsim_method must be 'errhmm' or 'qshmm' (got '${params.pbsim_method}')"
    }

    if (problems) {
        error "Parameter validation failed:\n  - " + problems.join('\n  - ')
        System.exit(1)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LONGREAD {
    validateFullRun()

    log.info """
    longread v${workflow.manifest.version}
  
    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

      Reference : ${params.bed} 
      Mapping   : ${params.transcript_gene}
      Sequence  : ${params.sequence}
      Model     : ${params.pbsim_model}
      Outdir    : ${params.outdir}
      Profile   : ${workflow.profile}
    """.stripIndent()

    ch_input = Channel.of(
        [ 
          [ id: params.prefix ], 
          file(params.bed, checkIfExists: true), 
          file(params.transcript_gene, checkIfExists: true) 
        ]
    )
    ch_sequence    = Channel.value(file(params.sequence, checkIfExists: true)).map { it -> [ it.baseName, it ] }
    ch_model       = Channel.value(file(params.pbsim_model, checkIfExists: true))

    SIMULATE(ch_input, ch_sequence, ch_model)
}

workflow {
    LONGREAD()
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nlongread finished. Results in: ${params.outdir}\n"
        : "\nlongread FAILED. See ${params.outdir}/pipeline_info for reports.\n")
}
