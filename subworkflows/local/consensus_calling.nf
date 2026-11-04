// Consensus across the callers that ran.
//
// SNVs: every selected caller contributes. Indels: only callers that report indels
// (MuTect1 is SNV-only). How callers are combined depends on how many take part (k):
//   k >= 3  n-1           a variant must be called by at least k-1 callers
//   k == 2  intersection  both callers must call it
//   k == 1  passthrough   the single caller's calls are returned as they are
// (k == 1 only happens for indels, e.g. MuTect1 + Strelka: SNVs are intersected and
// Strelka's indels pass through.) main.nf only runs this when at least two callers ran.

include { CONCAT_STRELKA    } from '../../modules/local/consensus'
include { ENSURE_INDEX      } from '../../modules/local/consensus'
include { SPLIT_SNVS_INDELS } from '../../modules/local/consensus'
include { CONSENSUS         } from '../../modules/local/consensus'
include { MERGE_CONSENSUS   } from '../../modules/local/consensus'

def consensusMethod(int k) {
    k >= 3 ? [name: 'n-1', min: k - 1]
           : k == 2 ? [name: 'intersection', min: 2]
                    : [name: 'passthrough', min: 1]
}

workflow CONSENSUS_CALLING {

    take:
    mutect_vcf        // (sampleName, vcf)            -- empty channel if MuTect1 did not run
    deepsomatic_vcf   // (sampleName, vcf)            -- empty channel if DeepSomatic did not run
    strelka_snvs      // (sampleName, vcf, tbi)       -- empty channel if Strelka did not run
    strelka_indels    // (sampleName, vcf, tbi)
    ref_fasta
    ref_fai

    main:

    // Callers that only report SNVs; they take no part in the indel consensus.
    def SNV_ONLY_CALLERS = ['mutect']

    // Callers that ran (take: values are channels, so read the list from params)
    def snv_callers   = params.callers.toString().tokenize(',').collect { it.trim() }.findAll { it }
    def indel_callers = snv_callers.findAll { !SNV_ONLY_CALLERS.contains(it) }

    def snv_plan   = consensusMethod(snv_callers.size())
    def indel_plan = indel_callers ? consensusMethod(indel_callers.size()) : null

    def snv_note = "SNV consensus: ${snv_plan.name} of ${snv_callers.join(', ')} (>= ${snv_plan.min} of ${snv_callers.size()} callers)"
    def indel_note = indel_plan
        ? "indel consensus: ${indel_plan.name} of ${indel_callers.join(', ')} (>= ${indel_plan.min} of ${indel_callers.size()} callers)"
        : "indel consensus: not run, ${snv_callers.join(', ')} only call SNVs"
    def notes = "callers: ${snv_callers.join(', ')} | ${snv_note} | ${indel_note}".toString()
    log.info notes.replace(' | ', '\n')

    // Strelka writes SNVs and indels as separate VCFs; join them so Strelka looks like every
    // other caller (one VCF per sample).
    strelka_ch = CONCAT_STRELKA(
        strelka_snvs
            .join(strelka_indels)
            .map { sample, snv_vcf, snv_tbi, indel_vcf, indel_tbi -> tuple(sample, snv_vcf, indel_vcf) }
    )

    // One (caller_id, vcf) per sample and caller, caller_id = "<sample>__<caller>" so tags and
    // file names stay unique.
    vcf_ch = mutect_vcf
        .map { sample, vcf -> tuple("${sample}__mutect".toString(), vcf) }
        .mix( deepsomatic_vcf.map { sample, vcf -> tuple("${sample}__deepsomatic".toString(), vcf) } )
        .mix( strelka_ch )

    indexed_ch = ENSURE_INDEX(vcf_ch)

    split_ch = SPLIT_SNVS_INDELS(indexed_ch, ref_fasta, ref_fai)

    // Group each sample's per-caller VCFs. A sample is released as soon as all of its callers
    // have finished.
    snv_in = split_ch.snvs
        .map { caller_id, vcf, tbi ->
            def sample_id = caller_id.split('__')[0]
            tuple(groupKey(sample_id, snv_callers.size()), vcf, tbi)
        }
        .groupTuple(by: 0)
        .map { sample_id, vcfs, tbis ->
            tuple(sample_id.toString(), 'snv', snv_plan.name, snv_plan.min, vcfs, tbis)
        }

    indel_in = Channel.empty()
    if (indel_plan) {
        indel_in = split_ch.indels
            .filter { caller_id, vcf, tbi -> indel_callers.contains(caller_id.split('__')[1]) }
            .map { caller_id, vcf, tbi ->
                tuple(groupKey(caller_id.split('__')[0], indel_callers.size()), vcf, tbi)
            }
            .groupTuple(by: 0)
            .map { sample_id, vcfs, tbis ->
                tuple(sample_id.toString(), 'indel', indel_plan.name, indel_plan.min, vcfs, tbis)
            }
    }

    consensus_ch = CONSENSUS(snv_in.mix(indel_in))

    // Per sample: the SNV consensus plus the indel consensus when there is one.
    n_parts = indel_plan ? 2 : 1
    merge_in = consensus_ch.consensus
        .map { sample_id, vcf, tbi, summary -> tuple(groupKey(sample_id, n_parts), vcf, tbi, summary) }
        .groupTuple(by: 0)
        .map { sample_id, vcfs, tbis, summaries -> tuple(sample_id.toString(), vcfs, tbis, summaries) }

    MERGE_CONSENSUS(merge_in, notes)

    emit:
    final_vcf = MERGE_CONSENSUS.out.merged            // (sampleName, vcf, tbi)
    consensus_log = MERGE_CONSENSUS.out.consensus_log
}
