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
    // reduce to top X loci by Fst
    //
    TOP_LOCI_FST(
        SNPIO_POPFILTER.out.filtered_vcf,
        SNPIO_POPFILTER.out.filtered_tbi,
        ch_joined_inputs.map { meta, vcf_file, tbi_file, popmap_file -> [meta, popmap_file] }
    )

    //
    // Join TOP_LOCI_FST output with popmap for proper matching
    //
    ch_top_vcf_with_popmap = TOP_LOCI_FST.out.top_vcf
        .join(ch_joined_inputs.map { meta, vcf_file, tbi_file, popmap_file -> [meta, popmap_file] })

    //
    // Simulate hybrid datasets (inspired by hybriddetective:simulate_gtfreq)
    //
    SIMULATE_HYBRIDS(
        ch_top_vcf_with_popmap.map { meta, vcf_file, popmap_file -> [meta, vcf_file] },
        ch_top_vcf_with_popmap.map { meta, vcf_file, popmap_file -> [meta, popmap_file] }
    )

    //
    // Prepare simulated data for input into NewHybrids
    //
    PREPARE_SIMULATION(
        ch_top_vcf_with_popmap.map { meta, vcf_file, popmap_file -> [meta, vcf_file] },
        ch_top_vcf_with_popmap.map { meta, vcf_file, popmap_file -> [meta, popmap_file] },
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
        ch_top_vcf_with_popmap.map { meta, vcf_file, popmap_file -> [meta, vcf_file] },
        ch_top_vcf_with_popmap.map { meta, vcf_file, popmap_file -> [meta, popmap_file] },
        Channel.value([[], []])
    )

    //
    // Run NewHybrids for the power analysis
    //
    RUN_NEWHYBRIDS(
        PREPARE_NEWHYBRIDS.out.newhybrids
    )
    ch_versions = ch_versions.mix( RUN_NEWHYBRIDS.out.versions )


    emit:
    versions     = ch_versions
    sim_trace    = POWER_ANALYSIS.out.pi_trace
    sim_map      = PREPARE_SIMULATION.out.index_map
    sim_result   = POWER_ANALYSIS.out.pofz
    nh_trace     = RUN_NEWHYBRIDS.out.pi_trace
    nh_map       = PREPARE_NEWHYBRIDS.out.index_map
    nh_result    = RUN_NEWHYBRIDS.out.pofz
}
