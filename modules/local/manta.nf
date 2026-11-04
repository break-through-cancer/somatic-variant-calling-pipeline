// Manta somatic SV / candidate-indel calling (tumor vs normal).
//
// Follows nf-core/modules manta/somatic. Strelka needs Manta's candidate small indels
// (--indelCandidates), so Manta always runs before Strelka. The SV VCFs are published too.
//
// Params used: outdir, wes.

// Convert (interval_list -> BED), compress and index the optional calling intervals. Manta and Strelka both require a
// bgzipped, tabix-indexed BED for --callRegions.
process BED_TABIX {
    tag "${bed.name}"
    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    path bed

    output:
    path "call_regions.bed.gz",     emit: bed_gz
    path "call_regions.bed.gz.tbi", emit: tbi

    script:
    // interval_list coordinates are 1-based; BED start is 0-based
    def shift = bed.name.toLowerCase() =~ /\.bed(\.gz)?$/ ? 0 : 1
    """
    set -euo pipefail
    case "${bed}" in
        *.gz) zcat ${bed} ;;
        *)    cat ${bed} ;;
    esac | grep -v -e '^#' -e '^@' | awk -F'\\t' -v OFS='\\t' '{ print \$1, \$2 - ${shift}, \$3 }' \\
         | sort -k1,1 -k2,2n | bgzip -c > call_regions.bed.gz
    tabix -p bed call_regions.bed.gz
    """
}


process MANTA_SOMATIC {
    tag "${sampleName}"
    label 'process_high'
    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
    publishDir "${params.outdir}/${sampleName}/manta", mode: 'copy'

    input:
    tuple val(sampleName), path(tumor), path(tumor_index), path(normal), path(normal_index)
    path fasta
    path fai
    path call_regions
    path call_regions_tbi

    output:
    tuple val(sampleName), path("${sampleName}.manta.candidate_small_indels.vcf.gz"), path("${sampleName}.manta.candidate_small_indels.vcf.gz.tbi"), emit: candidate_small_indels
    tuple val(sampleName), path("${sampleName}.manta.candidate_sv.vcf.gz"),          path("${sampleName}.manta.candidate_sv.vcf.gz.tbi"),          emit: candidate_sv
    tuple val(sampleName), path("${sampleName}.manta.diploid_sv.vcf.gz"),            path("${sampleName}.manta.diploid_sv.vcf.gz.tbi"),            emit: diploid_sv
    tuple val(sampleName), path("${sampleName}.manta.somatic_sv.vcf.gz"),            path("${sampleName}.manta.somatic_sv.vcf.gz.tbi"),            emit: somatic_sv

    script:
    // --exome only for targeted data; for WGS the call regions just restrict the contigs.
    def regions = call_regions ? "--callRegions ${call_regions}" : ''
    def exome   = params.wes.toString() == 'true' ? '--exome' : ''
    """
    set -euo pipefail

    configManta.py \\
        --tumorBam ${tumor} \\
        --normalBam ${normal} \\
        --reference ${fasta} \\
        ${regions} \\
        ${exome} \\
        --runDir manta

    python manta/runWorkflow.py -m local -j ${task.cpus}

    for t in candidateSmallIndels:candidate_small_indels candidateSV:candidate_sv diploidSV:diploid_sv somaticSV:somatic_sv; do
        src=\${t%%:*}
        dst=\${t##*:}
        mv manta/results/variants/\${src}.vcf.gz     ${sampleName}.manta.\${dst}.vcf.gz
        mv manta/results/variants/\${src}.vcf.gz.tbi ${sampleName}.manta.\${dst}.vcf.gz.tbi
    done
    """
}
