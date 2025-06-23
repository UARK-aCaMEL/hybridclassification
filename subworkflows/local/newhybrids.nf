//
// Run Steve Mussmann's Admixture Pipeline (AdmixPipe 3.0)
//

include { TABIX_BGZIP } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX } from '../../modules/nf-core/tabix/tabix/main'
include { SNPIO_POPFILTER } from '../../modules/local/snpio/pop_filter.nf'
include { TOP_LOCI_FST } from '../../modules/local/top_loci_fst.nf'
include { VCF_TO_GENEPOP } from '../../modules/local/vcf_to_genepop.nf'

workflow NEWHYBRIDS {
    take:
    vcf         // [ val(meta), *.vcf or *.vcf.gz ]
    tbi
    popmap

    main:
    ch_versions = Channel.empty()

    //
    // SNPio to remove loci not found in all groups
    //
    SNPIO_POPFILTER(
        vcf,
        tbi,
        popmap
    )

    //
    // reduce to top X loci by Fst
    //
    TOP_LOCI_FST(
        SNPIO_POPFILTER.out.filtered_vcf,
        SNPIO_POPFILTER.out.filtered_tbi,
        popmap
    )

    //
    // GENEPOP conversion
    //
    VCF_TO_GENEPOP(
        TOP_LOCI_FST.out.top_vcf,
        popmap
    )

    // Fetch results for the best K value
    // BESTK(
    //     CVSUM.out.cv_output,
    //     DISTRUCT.out.best_results
    // )
    // ch_versions = ch_versions.mix( BESTK.out.versions )

    emit:
    versions     = ch_versions
}
