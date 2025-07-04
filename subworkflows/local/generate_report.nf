include { PLOT_ADMIXTURE } from '../../modules/local/report/plot_admixture.nf'

workflow GENERATE_REPORT {
    take:
    report_inputs

    main:
    ch_versions = Channel.empty()
    ch_mqc_files = Channel.empty()

    //Admixture barplots
    PLOT_ADMIXTURE(
        report_inputs.map { m, i, p, k, bk -> tuple(m, k) },
        report_inputs.map { m, i, p, k, bk -> tuple(m, i) },
        report_inputs.map { m, i, p, k, bk -> tuple(m, p) }
    )
    ch_mqc_files = ch_mqc_files.mix( PLOT_ADMIXTURE.out.admixture_html )
    ch_versions = ch_versions.mix( PLOT_ADMIXTURE.out.versions )

    emit:
    mqc_files    = ch_mqc_files
    versions     = ch_versions
}
