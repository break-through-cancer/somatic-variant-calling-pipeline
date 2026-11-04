// MuTect1 scatter/gather with optional ContEst contamination estimation.
//
// This is the workflow body of the standalone mutect-1 pipeline, packaged as a subworkflow.
// MuTect1 needs BAM input (not CRAM) and one normal shared by all tumors.

include { INDEX_CONTEST_VCF } from '../../modules/local/mutect1'
include { CONTEST           } from '../../modules/local/mutect1'
include { SPLIT_INTERVALS   } from '../../modules/local/mutect1'
include { MUTECT1           } from '../../modules/local/mutect1'
include { GATHER_AND_FILTER } from '../../modules/local/mutect1'

workflow MUTECT1_CALLING {

    take:
    runs_ch      // (sampleName, tumor_bam, tumor_bai)
    normal_bam   // shared normal BAM
    normal_bai
    ref_fasta
    ref_fai
    ref_dict

    main:

    // Optional file inputs: an unset param is passed as [] so Nextflow stages nothing and the
    // process sees an empty string.
    def NO_FILE = []

    // force_calling_intervals only replaces target_list when force calling is explicitly
    // enabled; otherwise use target_list (or neither for genome-wide)
    use_force_calling_intervals = params.force_calling.toString() == 'true' && params.force_calling_intervals
    target_list = use_force_calling_intervals ? file(params.force_calling_intervals) :
                  (params.target_list ? file(params.target_list, checkIfExists: true) : NO_FILE)

    // Optional launch-form fields have been pre-filled with unrelated paths before (a BAM as
    // force_calling_intervals, an interval list as normal_panel). Fail up front with a clear
    // message rather than handing MuTect1 the wrong kind of file.
    expected_exts = [
        target_list      : ['.bed', '.interval_list', '.intervals'],
        dbsnp            : ['.vcf', '.vcf.gz'],
        cosmic           : ['.vcf', '.vcf.gz'],
        normal_panel     : ['.vcf', '.vcf.gz'],
        dbsnp_idx        : ['.idx', '.tbi'],
        cosmic_idx       : ['.idx', '.tbi'],
        normal_panel_idx : ['.idx', '.tbi']
    ]
    if (use_force_calling_intervals) {
        expected_exts.force_calling_intervals = ['.interval_list']
    }
    expected_exts.each { param_name, exts ->
        def value = params[param_name]
        if (value && !exts.any { value.toString().toLowerCase().endsWith(it) }) {
            error "Parameter ${param_name} must end in one of ${exts}, got: ${value}. " +
                  "Clear this field if you do not need it."
        }
    }

    // No fixed format for the read-group blacklist, so reject known-wrong file types
    wrong_blacklist_exts = ['.bam', '.bai', '.cram', '.crai', '.vcf', '.vcf.gz', '.interval_list',
                            '.bed', '.fasta', '.fa', '.fai', '.dict']
    if (params.read_group_blacklist &&
        wrong_blacklist_exts.any { params.read_group_blacklist.toString().toLowerCase().endsWith(it) }) {
        error "Parameter read_group_blacklist looks like the wrong file type, got: ${params.read_group_blacklist}. " +
              "Clear this field if you do not need it."
    }

    dbsnp    = params.dbsnp     ? file(params.dbsnp,     checkIfExists: true) : NO_FILE
    dbsnpIdx = params.dbsnp_idx ? file(params.dbsnp_idx, checkIfExists: true) : NO_FILE

    cosmic    = params.cosmic     ? file(params.cosmic,     checkIfExists: true) : NO_FILE
    cosmicIdx = params.cosmic_idx ? file(params.cosmic_idx, checkIfExists: true) : NO_FILE

    rgBlacklist = params.read_group_blacklist ? file(params.read_group_blacklist, checkIfExists: true) : NO_FILE

    normalPanel    = params.normal_panel     ? file(params.normal_panel,     checkIfExists: true) : NO_FILE
    normalPanelIdx = params.normal_panel_idx ? file(params.normal_panel_idx, checkIfExists: true) : NO_FILE

    // ContEst: estimate contamination per tumor, or use a fixed fraction
    if (params.contest_target_intervals && params.hapmap_vcf) {

        INDEX_CONTEST_VCF(file(params.hapmap_vcf, checkIfExists: true))

        snp6Bed = params.snp6_bed ? file(params.snp6_bed, checkIfExists: true) : NO_FILE

        CONTEST(
            runs_ch,
            normal_bam,
            normal_bai,
            file(params.contest_target_intervals, checkIfExists: true),
            ref_fasta,
            ref_fai,
            ref_dict,
            snp6Bed,
            INDEX_CONTEST_VCF.out.indexed_vcf
        )

        frac_ch = CONTEST.out.frac.map { pairName, fracFile ->
            tuple(pairName, fracFile.text.trim() as Float)
        }

    } else {

        frac_ch = runs_ch.map { pairName, t_bam, t_bai -> tuple(pairName, params.fracContam) }
    }

    runs_with_frac = runs_ch.join(frac_ch)

    // Scatter each tumor into interval shards
    SPLIT_INTERVALS(runs_ch, params.ref, target_list)

    shards_per_pair = SPLIT_INTERVALS.out.interval_files.transpose()

    mutect1_inputs = shards_per_pair
        .combine(runs_with_frac, by: 0)
        .map { pairName, intervalFile, t_bam, t_bai, fracContam ->
            tuple(pairName, intervalFile, fracContam, t_bam, t_bai)
        }

    MUTECT1(
        mutect1_inputs,
        normal_bam,
        normal_bai,
        ref_fasta,
        ref_fai,
        ref_dict,
        dbsnp,
        dbsnpIdx,
        cosmic,
        cosmicIdx,
        rgBlacklist,
        normalPanel,
        normalPanelIdx
    )

    // Gather the shards back into one call set per tumor
    gathered = MUTECT1.out.shard.groupTuple(by: [0, 1, 2])

    GATHER_AND_FILTER(gathered)

    emit:
    pass_vcf  = GATHER_AND_FILTER.out.pass_vcf_getzlab    // (sampleName, vcf)
    final_vcf = GATHER_AND_FILTER.out.final_vcf_neodisc   // (sampleName, vcf)
}
