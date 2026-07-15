/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
 * PROCESS_ISOSEQ: CCS per chunk -> merge CCS -> (optional) Iso-Seq clustering.
 */

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PBCCS }   from '../modules/local/pbccs/main.nf'
include { PBTK_PBMERGE as PBMERGE }         from '../modules/local/pbtk/pbmerge/main.nf'
include { ISOSEQ_CLUSTER2 } from '../modules/local/isoseq/cluster2/main.nf'
include { VALIDATE_CCS_CHUNKS } from '../modules/local/pacbio/validate_ccs_chunks/main.nf'
include { BAM_TO_FA } from '../modules/local/bamtofa/main.nf'

include { PUBLISH as PUBLISH_CCS } from '../modules/local/publish/main.nf'
include { PUBLISH as PUBLISH_CLUSTERED } from '../modules/local/publish/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROCESS_ISOSEQ {
    take:
      subreads            // [ meta, bam, pbi ] per chunk
      do_isoseq_cluster   // boolean

    main:
      ch_versions = Channel.empty()

      subreads
          .combine(Channel.of(1..params.ccs_chunk))   // INFO: cartesian product: N_bam × chunk combos
          .map { meta, bam, pbi, chunk_idx ->
              [ meta + [chunk: chunk_idx], bam, pbi ]
          }
          .set { ch_chunks }

      PBCCS(ch_chunks, params.ccs_chunk) // INFO: generate CCS from raw reads
      PBCCS.out.bam // INFO: update meta: update id (+chunkX) and store former id
      .map {
          def chunk   = it[0].chunk
          def parent  = it[0].id
          def child   = it[0].id + "." + chunk
          return [ [id:child, parent:parent, single_end:true, chunk:chunk], it[1] ]
      }
      .set { ch_pbccs_bam_updated }

      // INFO: group all chunks belonging to the same parent sample
      ch_pbccs_bam_updated
          .map { meta, bam -> [ meta.parent, meta, bam ] }   // INFO: key by parent
          .groupTuple(by: 0, size: params.ccs_chunk)              // INFO: wait for all chunks
          .map { parent, metas, bams ->
              // Reconstruct a clean meta for the merged output
              def meta_merged = [ id: parent, single_end: true ]
              [ meta_merged, bams ]                           // INFO: bams is now a List
          }
          .set { ch_pbccs_merged }

      // A CCS movie/ZMW may occur in exactly one chunk. Check this explicitly
      // so malformed input metadata fails with a useful error before pbmerge.
      VALIDATE_CCS_CHUNKS(ch_pbccs_merged)

      ch_pbccs_merged
          .join(VALIDATE_CCS_CHUNKS.out.valid, by: 0)
          .map { meta, bams, valid -> [ meta, bams ] }
          .set { ch_pbccs_validated }

      PBMERGE(ch_pbccs_validated) // INFO: merge chunks
      PUBLISH_CCS(PBMERGE.out.bam)

      ch_clustered = Channel.empty()
      ch_clustered_fa = Channel.empty()
      if (do_isoseq_cluster) {
          ISOSEQ_CLUSTER2(PBMERGE.out.bam)

          ch_clustered = ISOSEQ_CLUSTER2.out.bam
          PUBLISH_CLUSTERED(ISOSEQ_CLUSTER2.out.bam)

          BAM_TO_FA(ISOSEQ_CLUSTER2.out.bam)
          ch_clustered_fa = ch_clustered_fa.mix(BAM_TO_FA.out.singletons)
          ch_clustered_fa = ch_clustered_fa.mix(BAM_TO_FA.out.hq)

          ch_versions  = ch_versions.mix(ISOSEQ_CLUSTER2.out.versions)
      }

      ch_versions = ch_versions.mix(PBMERGE.out.versions)
      ch_versions = ch_versions.mix(PBCCS.out.versions)
      ch_versions = ch_versions.mix(VALIDATE_CCS_CHUNKS.out.versions)

    emit:
      ccs       = PBMERGE.out.bam
      clustered = ch_clustered
      fa        = ch_clustered_fa
      versions  = ch_versions
}
