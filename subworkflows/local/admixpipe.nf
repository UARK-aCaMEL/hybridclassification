//
// Run Steve Mussmann's Admixture Pipeline (AdmixPipe 3.0)
//

include { HTSLIB_BGZIPTABIX as DECOMPRESS_VCF } from '../../modules/nf-core/htslib/bgziptabix/main'
include { ADMIXTUREPIPELINE } from '../../modules/local/admixpipe/admixturepipeline.nf'
include { CLUMPAK }           from '../../modules/local/admixpipe/submitclumpak.nf'
include { CVSUM }             from '../../modules/local/admixpipe/cvsum.nf'
include { DISTRUCT }          from '../../modules/local/admixpipe/distructrerun.nf'
include { EVALADMIX }         from '../../modules/local/admixpipe/evaladmix.nf'
include { BESTK }             from '../../modules/local/bestK.nf'

workflow ADMIXPIPE {
    take:
    ch_input   // [ meta, vcf, popmap ]
    seed       // value: random seed for ADMIXTURE replicates

    main:
    ch_versions = Channel.empty()

    // Branch input VCF by extension
    ch_vcf_branch = ch_input.branch {
        vcfgz:     it[1].name.endsWith('.vcf.gz')
        vcf_plain: it[1].name.endsWith('.vcf')
    }

    // If input was vcf.gz, decompress
    DECOMPRESS_VCF(
        ch_vcf_branch.vcfgz
            .map { meta,vcf,popmap -> [ meta, vcf, [], [] ] },
        'decompress',
        false,
        'vcf'
    )

    // Combine uncompressed .vcf with decompressed .vcf
    ch_vcf = DECOMPRESS_VCF.out.output.mix(
        ch_vcf_branch.vcf_plain
            .map { meta,vcf,popmap -> tuple(meta,vcf) }
    )

    // Pass to ADMIXTURE pipeline
    ch_popmap = ch_input.map { meta,vcf,popmap -> tuple(meta,popmap) }
    ch_admixpipe_input = ch_vcf
        .join(ch_popmap)
        .map { meta,vcf,popmap -> tuple(meta,vcf,popmap) }
    ADMIXTUREPIPELINE( ch_admixpipe_input, seed )
    ch_versions = ch_versions.mix( ADMIXTUREPIPELINE.out.versions )

    // Run CLUMPAK
    clumpakIn  = ADMIXTUREPIPELINE.out.results
                    .join(ADMIXTUREPIPELINE.out.inds)
                    .join(ADMIXTUREPIPELINE.out.pops)
                    .join(ADMIXTUREPIPELINE.out.args_json)
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

    // Assess model fit with evalAdmix
    evalAdmixIn = ADMIXTUREPIPELINE.out.ped
                    .join(ADMIXTUREPIPELINE.out.map)
                    .join(ADMIXTUREPIPELINE.out.pfiles)
                    .join(ADMIXTUREPIPELINE.out.qfiles)
                    .join(ADMIXTUREPIPELINE.out.qfiles_json)
                    .join(ch_popmap)
                    .join(CLUMPAK.out.output)
                    .join(DISTRUCT.out.major_clusters)
                    .join(DISTRUCT.out.cvruns_json)
                    .join(DISTRUCT.out.qfilepaths_json)
    EVALADMIX( evalAdmixIn )
    ch_versions = ch_versions.mix( EVALADMIX.out.versions )

    // Fetch results for the best K value
    ch_bestk = CVSUM.out.cv_output
                .join( DISTRUCT.out.best_results )
    BESTK(
        ch_bestk.map { m, c, b -> [m, c] },
        ch_bestk.map { m, c, b -> [m, b] }
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
    qfilepaths   = DISTRUCT.out.qfilepaths_json
    corres       = EVALADMIX.out.corres
    fam          = EVALADMIX.out.fam
    versions     = ch_versions
}
