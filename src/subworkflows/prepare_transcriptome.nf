/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
 * PREPARE_TRANSCRIPTOME: prepare -> xloci extraction -> build PBSIM3 transcript input.
 */

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LONGREAD_PREPARE } from '../modules/local/longread/prepare/main.nf'
include { LONGREAD_PBSIM }   from '../modules/local/longread/pbsim3/main.nf'
include { XLOCI_EXON }    from '../modules/local/xloci/exon/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_TRANSCRIPTOME {
    take:
      ch_input        // [ meta, bed, transcript_gene ]
      ch_chrom_sizes  // path | []
      ch_sequence     // channel: [ meta, sequence ]

    main:
      ch_versions = Channel.empty()

      LONGREAD_PREPARE(ch_input, ch_chrom_sizes)
      ch_versions = ch_versions.mix(LONGREAD_PREPARE.out.versions)

      XLOCI_EXON(ch_sequence, LONGREAD_PREPARE.out.bed)
      ch_versions = ch_versions.mix(XLOCI_EXON.out.versions)

      // [ meta, sequences, isoform_depth, manifest ]
      ch_transcripts = Channel.empty()
      if (params.pbsim_mode == "trans") {
          ch_pbsim_in = XLOCI_EXON.out.tsv
              .join(LONGREAD_PREPARE.out.isoform_depth)
              .join(LONGREAD_PREPARE.out.manifest)

          LONGREAD_PBSIM(ch_pbsim_in)
          ch_transcripts = LONGREAD_PBSIM.out.transcript

          ch_versions = ch_versions.mix(LONGREAD_PBSIM.out.versions)
      } else if (params.pbsim_mode == "wgs") {
          ch_transcripts = XLOCI_EXON.out.fasta
      }

    emit:
      transcript    = ch_transcripts
      isoform_depth = LONGREAD_PREPARE.out.isoform_depth
      bed           = LONGREAD_PREPARE.out.bed
      stats         = LONGREAD_PREPARE.out.stats
      versions      = ch_versions
}
