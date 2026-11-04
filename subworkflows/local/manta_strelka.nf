// Manta -> Strelka2 somatic calling.
//
// Manta runs first and its candidate small indels are handed to Strelka, as in
// nf-core/sarek. An optional bgzipped+indexed BED restricts both callers to the regions
// given (used to drop decoy/alt contigs for WGS); pass [] for none.

include { MANTA_SOMATIC   } from '../../modules/local/manta'
include { STRELKA_SOMATIC } from '../../modules/local/strelka'

workflow MANTA_STRELKA {

    take:
    paired_ch        // (sampleName, tumor, tumor_index, normal, normal_index)
    normal_name      // normal sample name, used in the Strelka output names
    ref_fasta
    ref_fai
    regions          // bgzipped call-regions BED, or []
    regions_tbi      // its tabix index, or []

    main:

    MANTA_SOMATIC(paired_ch, ref_fasta, ref_fai, regions, regions_tbi)

    strelka_in = paired_ch.join(MANTA_SOMATIC.out.candidate_small_indels)

    STRELKA_SOMATIC(strelka_in, normal_name, ref_fasta, ref_fai, regions, regions_tbi)

    emit:
    snvs       = STRELKA_SOMATIC.out.snvs       // (sampleName, vcf, tbi)
    indels     = STRELKA_SOMATIC.out.indels     // (sampleName, vcf, tbi)
    somatic_sv = MANTA_SOMATIC.out.somatic_sv   // (sampleName, vcf, tbi)
}
