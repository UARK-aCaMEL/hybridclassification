/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ADMIXPIPE              } from '../subworkflows/local/admixpipe.nf'
include { NEWHYBRIDS             } from '../subworkflows/local/newhybrids.nf'
include { GENOMIC_CLINES         } from '../subworkflows/local/genomic_clines.nf'
include { GENERATE_REPORT        } from '../subworkflows/local/generate_report.nf'
include { SNPIO_FILTER           } from '../modules/local/snpio/pre_filter.nf'
include { SNPIO_SELECT           } from '../modules/local/snpio/select_pops.nf'
include { FIND_CANDIDATES        } from '../modules/local/find_candidates.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HYBRIDCLASSIFICATION {

    take:
    ch_vcf     // [meta, vcf]
    ch_tbi     // [meta, tbi]
    ch_popmap  // [meta, popmap]
    ch_speciesmap // [meta, speciesmap]
    ch_site_coords
    ch_geo_data
    ch_geo_data_dir
    ch_combinations

    main:

    ch_versions = Channel.empty()

    //
    // Generate subset VCFs for each test combination
    //
    ch_combinations
        .combine(ch_vcf)
        .combine(ch_tbi)
        .combine(ch_speciesmap)
        .map { pops, meta_vcf, vcf, meta_tbi, tbi, meta_popmap, popmap ->
            [pops, meta_vcf, vcf, tbi, popmap]
        }
        .set { ch_snpio_input }

    SNPIO_SELECT(
        ch_snpio_input.map { c, m, v, t, p -> [m, v] },
        ch_snpio_input.map { c, m, v, t, p -> [m, t] },
        ch_snpio_input.map { c, m, v, t, p -> [m, p] },
        ch_snpio_input.map { c, m, v, t, p -> c },
    )

    //
    //VCF pre-processing
    //
    SNPIO_SELECT.out.filtered_vcf
        .combine(ch_speciesmap)
        .map { pops, vcf, meta_popmap, popmap ->
            [pops, vcf, popmap]
        }
        .set { ch_filter_input }

    SNPIO_FILTER(
        ch_filter_input
    )
    ch_versions = ch_versions.mix(SNPIO_FILTER.out.versions)


    //
    // Run admixture pipeline on each filtered dataset
    //
    SNPIO_FILTER.out.filtered_vcf
        .combine( ch_speciesmap )
        .map { meta, vcf, popmap_meta, popmap -> tuple(meta, vcf, popmap) }
        .set { ch_admixpipe_input }

    ADMIXPIPE(
        ch_admixpipe_input
    )
    ch_versions = ch_versions.mix(ADMIXPIPE.out.versions)


    //
    // Generate inputs for newhybrids subworkflow
    //
    ch_joined  = ADMIXPIPE.out.k2_clumpp.join(ADMIXPIPE.out.inds)
    ch_find_in = ch_joined
        .combine(ch_popmap.map{meta, popmap -> popmap})
        .combine(ch_speciesmap.map{meta, popmap -> popmap})
    FIND_CANDIDATES(
        ch_find_in.map { m, k2, i, p, s -> tuple(m, k2) },
        ch_find_in.map { m, k2, i, p, s -> tuple(m, i) },
        ch_find_in.map { m, k2, i, p, s -> tuple(m, p) },
        ch_find_in.map { m, k2, i, p, s -> tuple(m, s) }
    )
    ch_versions = ch_versions.mix( FIND_CANDIDATES.out.versions )

    //
    // NewHybrids subworkflow
    //
    ch_joined_nh = SNPIO_FILTER.out.filtered_vcf
                        .join(SNPIO_FILTER.out.filtered_tbi)
                        .join(FIND_CANDIDATES.out.popmap)
    NEWHYBRIDS(
        ch_joined_nh.map { m, v, t, p -> tuple(m, v) },
        ch_joined_nh.map { m, v, t, p -> tuple(m, t) },
        ch_joined_nh.map { m, v, t, p -> tuple(m, p) }
    )
    ch_versions = ch_versions.mix( NEWHYBRIDS.out.versions )


    //
    // Optional: Genomic clines subworkflow
    //
    ch_bgc_text  = Channel.empty()
    if (params.run_bgc){
        ch_joined_bgc = SNPIO_FILTER.out.filtered_vcf
                        .join(SNPIO_FILTER.out.filtered_tbi)
                        .join(FIND_CANDIDATES.out.popmap)
        GENOMIC_CLINES(
            ch_joined_bgc.map { m, v, t, p -> tuple(m, v) },
            ch_joined_bgc.map { m, v, t, p -> tuple(m, t) },
            ch_joined_bgc.map { m, v, t, p -> tuple(m, p) }
        )
        ch_versions = ch_versions.mix( GENOMIC_CLINES.out.versions )
        ch_bgc_text  = GENOMIC_CLINES.out.bgc_text
    }

    //
    // Generate reports
    //

    GENERATE_REPORT(
        ADMIXPIPE.out.inds,
        ADMIXPIPE.out.pops,
        ADMIXPIPE.out.k2_clumpp,
        ADMIXPIPE.out.bestK_clumpp,
        NEWHYBRIDS.out.nh_result,
        NEWHYBRIDS.out.nh_trace,
        NEWHYBRIDS.out.nh_map,
        NEWHYBRIDS.out.sim_result,
        NEWHYBRIDS.out.sim_trace,
        NEWHYBRIDS.out.sim_map,
        NEWHYBRIDS.out.triangle_popmap,
        NEWHYBRIDS.out.triangle_hindex,
        NEWHYBRIDS.out.triangle_hindex_fixed,
        NEWHYBRIDS.out.masked_samples,
        FIND_CANDIDATES.out.popmap,
        ch_popmap,
        ch_speciesmap,
        ch_site_coords,
        ch_geo_data,
        ch_geo_data_dir,
        ch_bgc_text,
        ch_versions
    )


    emit:
    multiqc_report = GENERATE_REPORT.out.multiqc_report
    versions       = GENERATE_REPORT.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
