//
// Run Steve Mussmann's Admixture Pipeline (AdmixPipe 3.0)
//

include { TABIX_BGZIP }       from '../../modules/nf-core/tabix/bgzip/main'
include { ADMIXTUREPIPELINE } from '../../modules/local/admixpipe/admixturepipeline.nf'
include { CLUMPAK }           from '../../modules/local/admixpipe/submitclumpak.nf'
include { CVSUM }             from '../../modules/local/admixpipe/cvsum.nf'
include { DISTRUCT }          from '../../modules/local/admixpipe/distructrerun.nf'
include { BESTK }             from '../../modules/local/bestK.nf'

workflow ADMIXPIPE {
    take:
    ch_input   // [ meta, vcf, popmap ]

    main:
    ch_versions = Channel.empty()

    // Branch input VCF by extension
    ch_vcf_branch = ch_input.branch {
        vcfgz:     it[1].name.endsWith('.vcf.gz')
        vcf_plain: it[1].name.endsWith('.vcf')
    }

    // If input was vcf.gz, decompress
    TABIX_BGZIP(
        ch_vcf_branch.vcfgz
            .map { meta,vcf,popmap -> tuple(meta,vcf) }
    )
    ch_versions = ch_versions.mix( TABIX_BGZIP.out.versions )

    // Combine uncompressed .vcf with bgzipped .vcf
    ch_vcf = TABIX_BGZIP.out.output.mix(
        ch_vcf_branch.vcf_plain
            .map { meta,vcf,popmap -> tuple(meta,vcf) }
    )

    // Pass to ADMIXTURE pipeline
    ch_popmap = ch_input.map { meta,vcf,popmap -> tuple(meta,popmap) }
    ch_admixpipe_input = ch_vcf
        .join(ch_popmap)
        .map { meta,vcf,popmap -> tuple(meta,vcf,popmap) }
    ADMIXTUREPIPELINE( ch_admixpipe_input )
    ch_versions = ch_versions.mix( ADMIXTUREPIPELINE.out.versions )

    // Run CLUMPAK
    clumpakIn  = ADMIXTUREPIPELINE.out.results
                    .join(ADMIXTUREPIPELINE.out.inds)
                    .join(ADMIXTUREPIPELINE.out.pops)
    CLUMPAK( clumpakIn )
    ch_versions = ch_versions.mix( CLUMPAK.out.versions )

    // Run DISTRUCT
    distructIn  = ADMIXTUREPIPELINE.out.pfiles
                    .join(ADMIXTUREPIPELINE.out.qfiles)
                    .join(ADMIXTUREPIPELINE.out.pops)
                    .join(ADMIXTUREPIPELINE.out.inds)
                    .join(ADMIXTUREPIPELINE.out.logs)
                    .join(CLUMPAK.out.output)
    DISTRUCT( distructIn )
    ch_versions = ch_versions.mix( DISTRUCT.out.versions )

    // Compute best K from crossval
    cvsumIn  = DISTRUCT.out.cv.join(DISTRUCT.out.loglik)
    CVSUM( cvsumIn )
    ch_versions = ch_versions.mix( CVSUM.out.versions )

    // Fetch results for the best K value
    BESTK(
        CVSUM.out.cv_output,
        DISTRUCT.out.best_results
    )
    ch_versions = ch_versions.mix( BESTK.out.versions )

    emit:
    best_results = DISTRUCT.out.best_results
    bestK        = BESTK.out.bestK_file
    bestK_clumpp = BESTK.out.bestK_clumpp
    k2_clumpp    = BESTK.out.k2_clumpp
    inds         = ADMIXTUREPIPELINE.out.inds
    pops         = ADMIXTUREPIPELINE.out.pops
    cv_file      = CVSUM.out.cv_output
    versions     = ch_versions
}
