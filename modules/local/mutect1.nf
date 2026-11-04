// MuTect1 processes (scatter/gather + ContEst), taken from the standalone mutect-1 pipeline.
//
// Differences from the standalone version:
//   * GATHER_AND_FILTER emits (pairName, vcf) tuples so the VCFs can be joined to a sample
//     downstream (consensus).
//   * The force-calling-interval helper processes and the commented-out dedup processes are
//     not included.
//
// Params used here: outdir, scatter_count, padding_bp, downsample, max_mismatch_baseq_sum,
// force_calling, exclude_chimeric (defaults are set in nextflow.config).

process INDEX_CONTEST_VCF {
    tag "ContEst germline resource"
    label 'process_low'
    container "ghcr.io/jchen1095/contest:1.0.0"
    errorStrategy 'retry'
    maxRetries 2

    input:
    path hapmapVcf

    output:
    tuple path(hapmapVcf), path("${hapmapVcf}.tbi"), emit: indexed_vcf

    shell:
    '''
    set -euxo pipefail

    command -v tabix >/dev/null 2>&1 || {
        echo "ERROR: tabix is not installed in the ContEst container" >&2
        exit 1
    }

    tabix -f -p vcf !{hapmapVcf}
    test -s !{hapmapVcf}.tbi
    '''
}


// Port of getzlab's ContEst task, run once per tumor sample against the
// shared normal.
process CONTEST {
    tag "${pairName}"
    label 'process_high'
    container "ghcr.io/jchen1095/contest:1.0.0"
    publishDir "${params.outdir}/${pairName}/contest", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(pairName), path(t_bam), path(t_bai)
    path normal_bam
    path normal_bai
    path targetIntervals, name: 'contest_targets.interval_list'
    path ref_fasta
    path ref_fai
    path ref_dict
    path snp6Bed, name: 'snp6_targets.interval_list'
    tuple path(hapmapVcf), path(hapmapVcfIndex)

    output:
    tuple val(pairName), path('fraction_contamination.txt'), emit: frac
    path 'contamination.af.txt',                             emit: raw_af
    path 'contamination.base_report.txt',                    emit: base_report
    path 'contamination_validation.array_free.txt',           emit: validation

    shell:
    '''
    set -euxo pipefail

    if [ ! -f !{t_bam}.bai ]; then
        ln -s !{t_bai} !{t_bam}.bai
    fi

    if [ ! -f !{normal_bam}.bai ]; then
        ln -s !{normal_bai} !{normal_bam}.bai
    fi

    java_mem_mb=!{task.memory.toMega() - 1024}

    snp6_args=""
    if [ -n "!{snp6Bed}" ]; then
        snp6_args="-L !{snp6Bed} -isr INTERSECTION"
    fi

    /usr/local/jre1.7.0_71/bin/java \
        -Xmx${java_mem_mb}m \
        -Djava.io.tmpdir=$PWD \
        -jar /usr/local/bin/GenomeAnalysisTK.jar \
        -T ContEst \
        -I:eval !{t_bam} \
        -I:genotype !{normal_bam} \
        -L !{targetIntervals} \
        $snp6_args \
        -R !{ref_fasta} \
        -l INFO \
        -pf !{hapmapVcf} \
        -o contamination.af.txt \
        --trim_fraction 0.03 \
        --beta_threshold 0.05 \
        -br contamination.base_report.txt \
        -mbc 100 \
        --min_genotype_depth 30 \
        --min_genotype_ratio 0.8

    python /usr/local/bin/extract_contamination.py \
        contamination.af.txt \
        fraction_contamination.txt \
        contamination_validation.array_free.txt \
        !{pairName}

    test -s fraction_contamination.txt
    test -s contamination.af.txt
    test -s contamination.base_report.txt
    '''
}


// ---------------------------------------------------------------------------
// SPLIT INTERVALS
// ---------------------------------------------------------------------------
process SPLIT_INTERVALS {
    tag "${pairName}"
    label 'process_low'
    container "ghcr.io/jchen1095/split_intervals:v48-padded"
    publishDir "${params.outdir}/${pairName}/intervals", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 4

    input:
    tuple val(pairName), path(t_bam), path(t_bai)
    val ref
    path target_list

    output:
    tuple val(pairName), path('picard/*.interval_list'), emit: interval_files

    shell:
    '''
    set -euxo pipefail

    # Hardcoded for Ensembl/numeric style interval lists (Twist_HCEP_V1)

    selected_chrs=$(printf "chr%s," {1..22} X Y | sed 's/,$//')

    target_arg=""

    if [ -n "!{target_list}" ]; then
        target_arg="-target_list !{target_list}"
    fi

    python3 /app/split_intervals.py \
        -bam !{t_bam} \
        -bai !{t_bai} \
        -interval_type picard \
        -N !{params.scatter_count} \
        -padding !{params.padding_bp} \
        $target_arg \
        -chrs $selected_chrs

    for interval in picard/*.interval_list; do
        test -s "$interval"
    done
    '''
}

// ---------------------------------------------------------------------------
// MUTECT1
// ---------------------------------------------------------------------------
process MUTECT1 {
    tag "${pairName}:${task.index}"
    label 'process_medium'
    container "ghcr.io/jchen1095/mutect_1_getzlab:v36"
    stageInMode 'symlink'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(pairName), path(intervals), val(fracContam), path(t_bam), path(t_bai)

    path normal_bam
    path normal_bai

    path ref_fasta
    path ref_fai
    path ref_dict

    path dbsnp
    path dbsnpIdx

    path cosmic
    path cosmicIdx

    path readGroupBlacklist
    path normalPanel
    path normalPanelIdx

    output:
    tuple val(pairName),
          env('tumor_sample'),
          env('normal_sample'),
          path("${pairName}.MuTect1.call_stats.${task.index}.txt"),
          path("${pairName}.Mutect1.${task.index}.vcf"),
          path("${pairName}.MuTect1.shard.${task.index}.log"),
          emit: shard

    shell:
    '''
    set -euxo pipefail

    if [ ! -f !{t_bam}.bai ]; then
        ln -s !{t_bai} !{t_bam}.bai
    fi

    if [ ! -f !{normal_bam}.bai ]; then
        ln -s !{normal_bai} !{normal_bam}.bai
    fi

    tumor_sample=$(samtools view -H !{t_bam} \
        | awk -F'\t' '/^@RG/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SM:/) {
                    sub(/^SM:/, "", $i)
                    print $i
                }
            }
        }' \
        | sort -u)

    normal_sample=$(samtools view -H !{normal_bam} \
        | awk -F'\t' '/^@RG/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SM:/) {
                    sub(/^SM:/, "", $i)
                    print $i
                }
            }
        }' \
        | sort -u)

    [[ -n "$tumor_sample" ]] || {
        echo "ERROR: no SM tag in tumor BAM header" >&2
        exit 1
    }

    [[ -n "$normal_sample" ]] || {
        echo "ERROR: no SM tag in normal BAM header" >&2
        exit 1
    }

    [[ $(printf '%s\n' "$tumor_sample" | wc -l) -eq 1 ]] || {
        echo "ERROR: multiple tumor SM values: $tumor_sample" >&2
        exit 1
    }

    [[ $(printf '%s\n' "$normal_sample" | wc -l) -eq 1 ]] || {
        echo "ERROR: multiple normal SM values: $normal_sample" >&2
        exit 1
    }

    [[ "$tumor_sample" != "$normal_sample" ]] || {
        echo "ERROR: tumor and normal share the same SM value: $tumor_sample" >&2
        exit 1
    }

    java_mem_mb=!{task.memory.toMega() - 1024}

    extra_args=""

    if [ -n "!{dbsnp}" ]; then

        if [ -z "!{dbsnpIdx}" ]; then
            echo "ERROR: dbSNP VCF was supplied without dbsnp_idx" >&2
            exit 1
        fi

        if [ ! -e "!{dbsnp}.idx" ]; then
            ln -s "!{dbsnpIdx}" "!{dbsnp}.idx"
        fi

        test -e "!{dbsnp}.idx"

        extra_args="$extra_args --dbsnp !{dbsnp}"
    fi

    if [ -n "!{cosmic}" ]; then

        if [ -z "!{cosmicIdx}" ]; then
            echo "ERROR: COSMIC VCF was supplied without cosmic_idx" >&2
            exit 1
        fi

        if [ ! -e "!{cosmic}.idx" ]; then
            ln -s "!{cosmicIdx}" "!{cosmic}.idx"
        fi

        test -e "!{cosmic}.idx"

        extra_args="$extra_args --cosmic !{cosmic}"
    fi

    if [ -n "!{readGroupBlacklist}" ]; then
        extra_args="$extra_args --read_group_black_list !{readGroupBlacklist}"
    fi

    if [ -n "!{normalPanel}" ]; then
        extra_args="$extra_args --normal_panel !{normalPanel}"
    fi

    if [ "!{params.force_calling}" == "true" ]; then
        extra_args="$extra_args --force_output"
    fi

    if [ "!{params.exclude_chimeric}" == "true" ]; then
        extra_args="$extra_args --exclude_chimeric_reads"
    fi

    shard_log="!{pairName}.MuTect1.shard.!{task.index}.log"

    /usr/local/java/jdk1.7.0_80/bin/java \
        -Xmx${java_mem_mb}m \
        -jar /app/mutect.jar \
        --analysis_type MuTect \
        --tumor_sample_name "$tumor_sample" \
        -I:tumor !{t_bam} \
        -I:normal !{normal_bam} \
        --normal_sample_name "$normal_sample" \
        --reference_sequence !{ref_fasta} \
        --fraction_contamination !{fracContam} \
        --out !{pairName}.MuTect1.call_stats.!{task.index}.txt \
        --downsample_to_coverage !{params.downsample} \
        --max_read_mismatch_quality_score_sum !{params.max_mismatch_baseq_sum} \
        --vcf !{pairName}.Mutect1.!{task.index}.vcf \
        -L !{intervals} \
        $extra_args \
        > "$shard_log" 2>&1

    cat "$shard_log"

    if grep -Piq '(error)|(killed)|(java\\.lang\\.[a-zA-Z]*(exception|error):)' \
        "$shard_log"
    then
        echo "ERROR: MuTect1 shard !{task.index} for !{pairName} failed -- see log above" >&2
        exit 1
    fi

    test -s !{pairName}.MuTect1.call_stats.!{task.index}.txt
    test -s !{pairName}.Mutect1.!{task.index}.vcf
    test -s "$shard_log"
    '''
}


// ---------------------------------------------------------------------------
// GATHER + FILTER
// ---------------------------------------------------------------------------
process GATHER_AND_FILTER {
    tag "${pairName}"
    label 'process_medium'
    container "ghcr.io/jchen1095/mutect_1_getzlab:v36"

    // The PASS VCF is published as <pair>.mutect1.pass.vcf so downstream tools such as
    // consensus-calling can identify the caller by filename. saveAs only changes the
    // published copy; it is not part of the task cache key.
    publishDir "${params.outdir}/${pairName}", mode: 'copy',
        saveAs: { fn -> fn == "${pairName}.MuTect1.PASS.vcf" ? "${pairName}.mutect1.pass.vcf" : fn }

    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(pairName),
          val(tumorSampleName),
          val(normalSampleName),
          path(call_stats_files),
          path(vcf_files),
          path(log_files)

    output:
    path "${pairName}.MuTect1.call_stats.txt",         emit: call_stats
    path "${pairName}.PASS+TiNrisk.call_stats.txt",    emit: tin_call_stats
    path "${pairName}.MuTect1.vcf",                    emit: vcf
    tuple val(pairName), path("${pairName}.MuTect1.PASS.vcf"),    emit: pass_vcf_getzlab
    tuple val(pairName), path("${pairName}.mutectv1.final.vcf"), emit: final_vcf_neodisc
    path "${pairName}.MuTect1.log",                    emit: shard_logs

    shell:
    '''
    set -euxo pipefail

    printf '%s\n' !{call_stats_files} > call_stats.list
    printf '%s\n' !{vcf_files} > vcf.list
    printf '%s\n' !{log_files} > log.list

    merge_callstats.py \
        --FILE_ARRAY_PATH call_stats.list \
        !{pairName}.MuTect1.call_stats.txt

    merge_callstats.py \
        --FILE_ARRAY_PATH vcf.list \
        !{pairName}.MuTect1.orig.vcf

    : > !{pairName}.MuTect1.log
    while IFS= read -r logfile; do
        echo "===== $logfile =====" >> !{pairName}.MuTect1.log
        cat "$logfile" >> !{pairName}.MuTect1.log
    done < log.list

    python3 - <<PYCODE
vcf_file_in = "!{pairName}.MuTect1.orig.vcf"
vcf_file_out = "!{pairName}.MuTect1.vcf"

with open(vcf_file_in) as f_in, open(vcf_file_out, 'w') as f_out:
    for line in f_in:
        if line.startswith('#CHROM'):
            f_out.write('##normal_sample=!{normalSampleName}\\n')
            f_out.write('##tumor_sample=!{tumorSampleName}\\n')
        f_out.write(line)
PYCODE

    cat \
        <(sed -n '/^#/p' !{pairName}.MuTect1.vcf) \
        <(
            awk -F '\t' \
                '$7 == "PASS" && $1 ~ /[0-9]+/' \
                !{pairName}.MuTect1.vcf \
                | sort -k1,1n -k2,2n
        ) \
        <(
            awk -F '\t' \
                '$7 == "PASS" && $1 !~ /[0-9]+/' \
                !{pairName}.MuTect1.vcf \
                | sort -k1,1 -k2,2n
        ) \
        > !{pairName}.MuTect1.PASS.vcf

    cat \
        <(sed -n '/^#/p' !{pairName}.MuTect1.vcf) \
        <(
            awk -F '\t' \
                '!/^#/ && !/REJECT/ && $1 ~ /[0-9]+/' \
                !{pairName}.MuTect1.vcf \
                | sort -k1,1n -k2,2n
        ) \
        <(
            awk -F '\t' \
                '!/^#/ && !/REJECT/ && $1 !~ /[0-9]+/' \
                !{pairName}.MuTect1.vcf \
                | sort -k1,1 -k2,2n
        ) \
        > !{pairName}.mutectv1.final.vcf

    awk -F '\t' \
        'NR <= 2 ||
         $51 == "KEEP" ||
         $50 == "alt_allele_in_normal" ||
         $50 == "normal_lod" ||
         $50 == "normal_lod,alt_allele_in_normal" ||
         $50 == "alt_allele_in_normal,normal_lod"' \
        !{pairName}.MuTect1.call_stats.txt \
        > !{pairName}.PASS+TiNrisk.call_stats.txt

    test -s !{pairName}.MuTect1.call_stats.txt
    test -s !{pairName}.MuTect1.vcf
    test -s !{pairName}.MuTect1.PASS.vcf
    test -s !{pairName}.mutectv1.final.vcf
    test -s !{pairName}.MuTect1.log
    '''
}
