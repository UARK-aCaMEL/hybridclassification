//
// Run Steve Mussmann's Admixture Pipeline (AdmixPipe 3.0)
//

include { TABIX_BGZIP } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX } from '../../modules/nf-core/tabix/tabix/main'
include { SNPIO_POPFILTER } from '../../modules/local/snpio/pop_filter.nf'
include { TOP_LOCI_FST } from '../../modules/local/top_loci_fst.nf'
include { VCF_TO_NEWHYBRIDS as PREPARE_SIMULATION } from '../../modules/local/vcf_to_newhybrids.nf'
include { VCF_TO_NEWHYBRIDS as PREPARE_NEWHYBRIDS } from '../../modules/local/vcf_to_newhybrids.nf'
include { SIMULATE_HYBRIDS } from '../../modules/local/hybrid_simulator.nf'
include { RUN_NEWHYBRIDS as POWER_ANALYSIS } from '../../modules/local/newhybrids.nf'
include { RUN_NEWHYBRIDS } from '../../modules/local/newhybrids.nf'

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
    // Simulate hybrid datasets (inspired by hybriddetective:simulate_gtfreq)
    //
    SIMULATE_HYBRIDS(
        TOP_LOCI_FST.out.top_vcf,
        popmap
    )

    //
    // Prepare simulated data for input into NewHybrids
    //
    PREPARE_SIMULATION(
        TOP_LOCI_FST.out.top_vcf,
        popmap,
        SIMULATE_HYBRIDS.out.vcf
    )

    //
    // Run NewHybrids for the power analysis
    //
    POWER_ANALYSIS(
        PREPARE_SIMULATION.out.newhybrids
    )
    ch_versions = ch_versions.mix( POWER_ANALYSIS.out.versions )

    //
    // Prepare full dataset for NewHybrids analysis
    //
    PREPARE_NEWHYBRIDS(
        TOP_LOCI_FST.out.top_vcf,
        popmap,
        [[], []]
    )

    //
    // Run NewHybrids for the power analysis
    //
    RUN_NEWHYBRIDS(
        PREPARE_NEWHYBRIDS.out.newhybrids
    )
    ch_versions = ch_versions.mix( POWER_ANALYSIS.out.versions )


    emit:
    versions     = ch_versions
    // sim_trace
    // sim_map
    // sim_result
    // nh_trace
    // nh_map
    // nh_result
}
