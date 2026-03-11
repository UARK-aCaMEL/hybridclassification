process BGC_HM {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/tkchafin/bgchm-r:1.0'

    input:
    tuple val(meta),  path(glikp0)
    tuple val(meta2), path(glikp1)
    tuple val(meta3), path(glikadmix)
    tuple val(meta4), path(glikp0_samples)
    tuple val(meta5), path(glikp1_samples)
    tuple val(meta6), path(glikadmix_samples)
    tuple val(meta7), path(locus_order)

    output:
    path "versions.yml",                                      emit: versions
    tuple val(meta), path("results_${meta.id}/plots"),        emit: bgc_plots
    tuple val(meta), path("results_${meta.id}/text"),         emit: bgc_text
    tuple val(meta), path("results_${meta.id}/rdata"),        emit: bgc_rdata

    script:
    def args = task.ext.args ?: ''
    """
    run_bgc_pipeline.R \\
        --p0 ${glikp0} \\
        --p1 ${glikp1} \\
        --admix ${glikadmix} \\
        --outdir results_${meta.id} \\
        --prefix ${meta.id} \\
        --n_iters ${params.bgc_iters} \\
        --p_warmup ${params.bgc_burnin} \\
        --ci 0.90 \\
        --n_thin ${params.bgc_thin} \\
        ${args}

    # Copy the order files into results_*/text using the same naming style as R outputs
    cp ${glikp0_samples}    results_${meta.id}/text/${meta.id}__input_order__GlikP0_samples.txt
    cp ${glikp1_samples}    results_${meta.id}/text/${meta.id}__input_order__GlikP1_samples.txt
    cp ${glikadmix_samples} results_${meta.id}/text/${meta.id}__input_order__GlikADMIX_samples.txt
    cp ${locus_order}       results_${meta.id}/text/${meta.id}__input_order__locus_order.txt

    # Log versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1 | sed -E 's/^R version ([0-9.]+).*/\\1/')
        bgchm: \$(Rscript -e 'pkgs <- c("bgchm","bgc_hm","bgc-hm","bgc.hm","bgc","bgcTools"); v <- NA; for(p in pkgs){ if (requireNamespace(p, quietly=TRUE)) { v <- as.character(utils::packageVersion(p)); break } }; if (is.na(v)) v <- "NA"; cat(v)')
    END_VERSIONS
    """
}
