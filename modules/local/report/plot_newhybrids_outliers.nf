process NH_PLOT_OUTLIERS {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
    tuple val(meta), path(hindex)
    tuple val(meta2), path(hindex_popmap)
    tuple val(meta3), path(pofz)
    tuple val(meta3), path(index_map)

    output:
    tuple val(meta), path("${meta.id}_nh_outliers_mqc.html"),     emit: plot_html

    script:
    def args = task.ext.args ?: ''
    """
    plot_nh_outliers.py \\
        --hindex ${hindex} \\
        --nh_results ${pofz} \\
        --nh_index ${index_map} \\
        --hindex_popmap ${hindex_popmap} \\
        --prob_threshold ${params.prob_threshold} \\
        --alpha ${params.outlier_alpha} \\
        --template ${baseDir}/assets/multiqc_nh_outliers.html \\
        --out "${meta.id}_nh_outliers_mqc.html" \\
        ${args}
    """
}