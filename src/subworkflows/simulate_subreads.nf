/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
 * SIMULATE_SUBREADS: split -> PBSIM3 per chunk -> merge into one subreads BAM.
 */

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LONGREAD_SPLIT } from '../modules/local/longread/split/main.nf'
include { PBSIM3 }         from '../modules/local/pbsim3/main.nf'
include { PBTK_PBMERGE as PBMERGE }        from '../modules/local/pbtk/pbmerge/main.nf'
include { NORMALIZE_PACBIO_RG }             from '../modules/local/pacbio/normalize_rg/main.nf'
include { PUBLISH as PUBLISH_MAF } from '../modules/local/publish/main.nf'
include { PUBLISH as PUBLISH_REF_MAP } from '../modules/local/publish/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SIMULATE_SUBREADS {
    take:
      ch_transcript  // [ meta, pbsim.transcript.tsv ] | [ meta, transcript.fasta ]
      ch_model       // path (PBSIM3 model)

    main:
      ch_versions = Channel.empty()

      if (params.pbsim_mode == "trans") {
          LONGREAD_SPLIT(ch_transcript)
          ch_versions = ch_versions.mix(LONGREAD_SPLIT.out.versions)

          // Fan out to one [ meta, chunk_tsv ] per chunk (normalise single-match to a list first).
          ch_chunks = LONGREAD_SPLIT.out.chunks
              .map { meta, files -> [ meta, files instanceof List ? files : [ files ] ] }
              .transpose()
      } else if (params.pbsim_mode == "wgs") {
          ch_transcript
              .splitFasta(
                  by: params.pbsim_chunks,
                  file: 'chunk'
              )
              .map { meta, chunk ->
                  def chunk_id = (chunk.name =~ /(chunk\.\d+)/)[0][1]
                  def updated_meta = meta + [
                      id       : "${meta.id}.${chunk_id}",
                      id_former: meta.id
                  ]
                  [ updated_meta, chunk ]
              }
              .set { ch_chunks }
      }

      PBSIM3(ch_chunks, ch_model)
      ch_versions = ch_versions.mix(PBSIM3.out.versions)

      // Collect all raw subread BAMs. The normalizer presents them as one
      // synthetic movie with globally unique ZMWs, as required by CCS chunking.
      PBSIM3.out.bam
        .map { meta, bams ->
            def parent = meta.id_former ?: meta.id
            [ parent, bams instanceof List ? bams : [ bams ] ]
        }
        .groupTuple(by: 0)
        .map { parent, bam_groups -> [ [ id: parent ], bam_groups.flatten() ] }
        .set { ch_bam_raw_grouped }

      NORMALIZE_PACBIO_RG(ch_bam_raw_grouped)
      ch_versions = ch_versions.mix(NORMALIZE_PACBIO_RG.out.versions)

      NORMALIZE_PACBIO_RG.out.bam
        .map { meta, bams -> [ meta, bams instanceof List ? bams : [bams] ] }
        .set { ch_bam_grouped }

      // pbmerge preserves PacBio sort order and creates the PBI consumed by CCS.
      PBMERGE(ch_bam_grouped)
      ch_versions = ch_versions.mix(PBMERGE.out.versions)

      PBMERGE.out.pbi
        .join(
            PBMERGE.out.bam,
            by: 0   // join on meta
        )
        .map { meta, pbi, bam ->
            def meta_updated = meta + [ indexed: true ]
            [ meta_updated, bam, pbi ]
        }.set { ch_bam }

     PBSIM3.out.maf
        .map { meta, maf -> maf }
        .collect()
        .map { mafs -> [ [ id : 'merged' ], mafs] }
        .set { ch_mafs_grouped }
     PUBLISH_MAF(ch_mafs_grouped)

     PBSIM3.out.map
        .map { meta, map -> map }
        .collect()
        .map { maps -> [ [ id : 'merged' ], maps] }
        .set { ch_maps_grouped }
     PUBLISH_REF_MAP(ch_maps_grouped)

    emit:
      bam       = ch_bam  // [ meta, ${prefix}.subreads.bam, ${prefix}.subreads.pbi ]
      zmw_map   = NORMALIZE_PACBIO_RG.out.map
      versions  = ch_versions
}
