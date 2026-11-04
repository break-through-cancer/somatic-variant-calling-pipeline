#!/usr/bin/env nextflow
/*
 * Somatic variant calling with optional consensus.
 *
 * Input: preprocessed tumor/normal BAM or CRAM files (one shared normal). No alignment is done.
 * Callers (pick any): MuTect1, DeepSomatic, Manta -> Strelka2.
 * Consensus: when run_consensus is on and at least two callers were selected, the callers'
 * calls are combined into a consensus VCF on top of each caller's own results.
 *
 * params.runs, normal_sample_name, normal_bam, normal_bam_index and callers are set by
 * .cirro/preprocess.py from the input dataset; the rest come from the launch form.
 */
nextflow.enable.dsl = 2

include { MUTECT1_CALLING   } from './subworkflows/local/mutect1_calling'
include { DEEPSOMATIC       } from './modules/local/deepsomatic'
include { BED_TABIX         } from './modules/local/manta'
include { MANTA_STRELKA     } from './subworkflows/local/manta_strelka'
include { CONSENSUS_CALLING } from './subworkflows/local/consensus_calling'

workflow {

    def KNOWN_CALLERS = ['mutect', 'deepsomatic', 'strelka']

    // -----------------------------------------------------------------------
    // Validate inputs
    // -----------------------------------------------------------------------
    if (!params.runs) {
        error "params.runs is empty -- did the Cirro preprocess.py hook run? (see .cirro/preprocess.py)"
    }
    if (!params.normal_bam || !params.normal_bam_index) {
        error "params.normal_bam / normal_bam_index are missing -- a matched normal is required"
    }
    if (!params.ref_fasta || !params.ref_fai) {
        error "Missing required reference params: ref_fasta, ref_fai"
    }

    def callers = params.callers.toString().tokenize(',').collect { it.trim() }.findAll { it }
    if (!callers) {
        error "No callers selected; tick at least one of MuTect1, DeepSomatic, Strelka"
    }
    def unknown = callers - KNOWN_CALLERS
    if (unknown) {
        error "Unknown caller(s): ${unknown.join(', ')}. Valid: ${KNOWN_CALLERS.join(', ')}"
    }

    def run_consensus = params.run_consensus.toString() == 'true'
    if (run_consensus && callers.size() < 2) {
        log.warn "run_consensus is on but only one caller (${callers[0]}) was selected; skipping consensus"
        run_consensus = false
    }

    if ('mutect' in callers) {
        if (!params.ref_dict) {
            error "MuTect1 needs ref_dict"
        }
        // MuTect1 (GATK3) reads BAM only
        def cram = (params.runs.collect { it.tumor_bam } + [params.normal_bam]).findAll { it.toString().endsWith('.cram') }
        if (cram) {
            error "MuTect1 cannot read CRAM, but the input has CRAM files (${cram.take(3).join(', ')}). " +
                  "Untick MuTect1 or provide BAMs."
        }
    }

    log.info "callers: ${callers.join(', ')} | consensus: ${run_consensus ? 'yes' : 'no'} | tumor samples: ${params.runs.size()}"

    // -----------------------------------------------------------------------
    // Inputs
    // -----------------------------------------------------------------------
    ref_fasta = file(params.ref_fasta, checkIfExists: true)
    ref_fai   = file(params.ref_fai,   checkIfExists: true)
    ref_dict  = params.ref_dict ? file(params.ref_dict, checkIfExists: true) : []

    normal_bam = file(params.normal_bam,       checkIfExists: true)
    normal_bai = file(params.normal_bam_index, checkIfExists: true)

    // (sampleName, tumor, tumor_index, normal, normal_index), one row per tumor
    paired_ch = Channel.fromList(
        params.runs.collect { r ->
            tuple(
                r.sample_name,
                file(r.tumor_bam,       checkIfExists: true),
                file(r.tumor_bam_index, checkIfExists: true),
                normal_bam,
                normal_bai
            )
        }
    )

    // -----------------------------------------------------------------------
    // Callers
    // -----------------------------------------------------------------------
    mutect_vcf      = Channel.empty()
    deepsomatic_vcf = Channel.empty()
    strelka_snvs    = Channel.empty()
    strelka_indels  = Channel.empty()

    if ('mutect' in callers) {
        MUTECT1_CALLING(
            paired_ch.map { name, t_bam, t_bai, n_bam, n_bai -> tuple(name, t_bam, t_bai) },
            normal_bam,
            normal_bai,
            ref_fasta,
            ref_fai,
            ref_dict
        )
        mutect_vcf = MUTECT1_CALLING.out.pass_vcf
    }

    if ('deepsomatic' in callers) {
        DEEPSOMATIC(paired_ch, ref_fasta, ref_fai)
        deepsomatic_vcf = DEEPSOMATIC.out.vcf
    }

    if ('strelka' in callers) {
        // Optional call regions for Manta/Strelka (bgzipped + indexed first)
        if (params.call_regions) {
            BED_TABIX(file(params.call_regions, checkIfExists: true))
            regions     = BED_TABIX.out.bed_gz
            regions_tbi = BED_TABIX.out.tbi
        } else {
            regions     = []
            regions_tbi = []
        }

        MANTA_STRELKA(
            paired_ch,
            params.normal_sample_name,
            ref_fasta,
            ref_fai,
            regions,
            regions_tbi
        )
        strelka_snvs   = MANTA_STRELKA.out.snvs
        strelka_indels = MANTA_STRELKA.out.indels
    }

    // -----------------------------------------------------------------------
    // Consensus
    // -----------------------------------------------------------------------
    if (run_consensus) {
        CONSENSUS_CALLING(
            mutect_vcf,
            deepsomatic_vcf,
            strelka_snvs,
            strelka_indels,
            ref_fasta,
            ref_fai
        )
    }
}
