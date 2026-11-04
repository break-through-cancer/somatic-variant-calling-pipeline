// DeepSomatic (NVIDIA Parabricks) -- taken unchanged from the standalone deepsomatic pipeline.
// Params used: outdir, num_streams_per_gpu, enable_gvcf.

process DEEPSOMATIC {
    tag "${sampleName}"
    label 'process_high'
    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"
    // Published files carry the caller name (<sample>.deepsomatic.vcf.gz) so downstream
    // tools such as consensus-calling can tell callers apart by filename. saveAs only
    // changes the published copy; it is not part of the task cache key.
    publishDir "${params.outdir}/${sampleName}", mode: 'copy',
        saveAs: { fn -> fn.replace("${sampleName}.", "${sampleName}.deepsomatic.") }
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sampleName), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(sampleName), path("${sampleName}.vcf.gz"), emit: vcf
    path "${sampleName}.vcf.gz.tbi", optional: true, emit: tbi
    path "${sampleName}.log", emit: log

    shell:
    '''
    set -euxo pipefail

    # Run DeepSomatic
    pbrun deepsomatic \
        --in-tumor-bam !{tumor_bam} \
        --in-normal-bam !{normal_bam} \
        --ref !{ref_fasta} \
        --out-variants !{sampleName}.vcf.gz \
        --num-streams-per-gpu !{params.num_streams_per_gpu} \
        !{params.enable_gvcf ? '--run-deepsomatic-gvcf' : ''} \
        > !{sampleName}.log 2>&1

    # Verify output. pbrun may not write a .tbi, so index it if we can, but don't fail the
    # task over a missing index (the VCF itself can take hours to produce).
    test -s !{sampleName}.vcf.gz
    if [ ! -s !{sampleName}.vcf.gz.tbi ] && command -v tabix >/dev/null 2>&1; then
        tabix -p vcf !{sampleName}.vcf.gz || true
    fi
    '''
}
