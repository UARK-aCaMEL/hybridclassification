//
// Run Steve Mussmann's Admixture Pipeline (AdmixPipe 3.0)
//

include { TABIX_BGZIP } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX } from '../../modules/nf-core/tabix/tabix/main'
include { SNPIO_POPFILTER } from '../../modules/local/snpio/pop_filter.nf'
include { VCF2BGC } from '../../modules/local/vcf2bgc.nf'
include { BGC_HM } from '../../modules/local/run_bgchm.nf'

workflow GENOMIC_CLINES {
    take:
    vcf         // [ val(meta), *.vcf or *.vcf.gz ]
    tbi         // [ val(meta), *.tbi ]
    popmap      // [ val(meta), popmap ]

    main:
    ch_versions = Channel.empty()

    // Join vcf, tbi, and popmap by meta to ensure proper matching
    ch_joined_inputs = vcf.join(tbi).join(popmap)

    //
    // SNPio to remove loci not found in all groups
    //
    SNPIO_POPFILTER(
        ch_joined_inputs.map { meta, vcf_file, tbi_file, popmap_file -> [meta, vcf_file] },
        ch_joined_inputs.map { meta, vcf_file, tbi_file, popmap_file -> [meta, tbi_file] },
        ch_joined_inputs.map { meta, vcf_file, tbi_file, popmap_file -> [meta, popmap_file] }
    )

    //
    // Format conversion and filter on delta
    //
    ch_top_loci = ch_joined_inputs
        .map{ meta, vcf_file, tbi_file, popmap_file -> [meta, popmap_file] }
        .join( SNPIO_POPFILTER.out.filtered_vcf )
        .join( SNPIO_POPFILTER.out.filtered_tbi )
    VCF2BGC(
        ch_top_loci.map { m, p, v, t -> [m, v] },
        ch_top_loci.map { m, p, v, t -> [m, p] }
    )

    //
    // Run BGC
    //
    BGC_HM(
        VCF2BGC.out.parent1,
        VCF2BGC.out.parent2,
        VCF2BGC.out.admix,
        VCF2BGC.out.parent1_samples,
        VCF2BGC.out.parent2_samples,
        VCF2BGC.out.admix_samples,
        VCF2BGC.out.locus_order
    )
    ch_versions = ch_versions.mix( BGC_HM.out.versions )


    emit:
    versions  = ch_versions
    bgc_plots = BGC_HM.out.bgc_plots
    bgc_text  = BGC_HM.out.bgc_text
    bgc_rdata = BGC_HM.out.bgc_rdata

}
