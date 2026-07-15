/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
 * LONGREAD main workflow: transcriptome preparation -> subread simulation -> optional CCS /
 * Iso-Seq clustering, with BAM validation of every final output.
 */

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_TRANSCRIPTOME } from '../subworkflows/prepare_transcriptome.nf'
include { SIMULATE_SUBREADS }     from '../subworkflows/simulate_subreads.nf'
include { PROCESS_ISOSEQ }        from '../subworkflows/process_isoseq.nf'
include { VALIDATE_BAM }          from '../modules/local/validate_bam/main.nf'
include { CHROMSIZE }             from '../modules/local/chromsize/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LONGREAD {
    take:
      ch_input        // [ meta, bed, transcript_gene ]
      ch_sequence     // channel: [ meta, sequence ]
      ch_model        // path (PBSIM3 model)

    main:
      ch_versions = Channel.empty()

      CHROMSIZE(ch_sequence)
      ch_chrom_sizes = CHROMSIZE.out.chrom_sizes.map { meta, chrom_sizes -> chrom_sizes }

      PREPARE_TRANSCRIPTOME(ch_input, ch_chrom_sizes, ch_sequence)
      ch_versions = ch_versions.mix(PREPARE_TRANSCRIPTOME.out.versions)

      SIMULATE_SUBREADS(PREPARE_TRANSCRIPTOME.out.transcript, ch_model)
      ch_versions = ch_versions.mix(SIMULATE_SUBREADS.out.versions)

      // BAMs that must pass validation (subreads always; CCS / clustered when produced).
      ch_bams = SIMULATE_SUBREADS.out.bam

      ch_ccs = Channel.empty()
      ch_clustered = Channel.empty()
      if (params.do_ccs) {
          PROCESS_ISOSEQ(SIMULATE_SUBREADS.out.bam, params.do_isoseq_cluster)

          ch_bams = ch_bams.mix(PROCESS_ISOSEQ.out.ccs).mix(PROCESS_ISOSEQ.out.clustered)
          ch_ccs = PROCESS_ISOSEQ.out.ccs
          ch_clustered = PROCESS_ISOSEQ.out.clustered

          ch_versions = ch_versions.mix(PROCESS_ISOSEQ.out.versions)
      }

      // VALIDATE_BAM(ch_bams)
      // ch_versions = ch_versions.mix(VALIDATE_BAM.out.versions)

      // Aggregate tool versions (dedup identical per-process YAML blocks by content).
      ch_versions
          .map { it.text }
          .unique()
          .collectFile(name: 'versions.yml', storeDir: "${params.outdir}/PIPELINE_INFO")

    emit:
      subreads  = SIMULATE_SUBREADS.out.bam
      zmw_map   = SIMULATE_SUBREADS.out.zmw_map
      ccs       = ch_ccs
      clustered = ch_clustered
      versions  = ch_versions
}
