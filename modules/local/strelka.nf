// Strelka2 somatic SNV / indel calling (tumor vs normal), seeded with Manta's candidate
// small indels.
//
// Follows nf-core/modules strelka/somatic. Outputs are published as
// <tumor>_vs_<normal>.strelka.somatic_snvs.vcf.gz / .somatic_indels.vcf.gz, the same names
// nf-core/sarek produces, so results look the same as earlier sarek-based runs.
//
// Params used: outdir, wes.

process STRELKA_SOMATIC {
    tag "${sampleName}"
    label 'process_high'
    container 'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1'
    publishDir "${params.outdir}/${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path(tumor), path(tumor_index), path(normal), path(normal_index), path(manta_indels), path(manta_indels_tbi)
    val normal_name
    path fasta
    path fai
    path call_regions
    path call_regions_tbi

    output:
    tuple val(sampleName), path("${prefix}.strelka.somatic_snvs.vcf.gz"),   path("${prefix}.strelka.somatic_snvs.vcf.gz.tbi"),   emit: snvs
    tuple val(sampleName), path("${prefix}.strelka.somatic_indels.vcf.gz"), path("${prefix}.strelka.somatic_indels.vcf.gz.tbi"), emit: indels

    script:
    prefix = "${sampleName}_vs_${normal_name}"
    def regions = call_regions ? "--callRegions ${call_regions}" : ''
    def exome   = params.wes.toString() == 'true' ? '--exome' : ''
    """
    set -euo pipefail

    configureStrelkaSomaticWorkflow.py \\
        --tumorBam ${tumor} \\
        --normalBam ${normal} \\
        --referenceFasta ${fasta} \\
        --indelCandidates ${manta_indels} \\
        ${regions} \\
        ${exome} \\
        --runDir strelka

    python strelka/runWorkflow.py -m local -j ${task.cpus}

    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.strelka.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.strelka.somatic_snvs.vcf.gz.tbi
    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.strelka.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.strelka.somatic_indels.vcf.gz.tbi
    """
}
