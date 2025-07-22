process PLOT_TRIANGLE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(hindex)
        tuple val(meta2), path(hindex_fixed)
        tuple val(meta3), path(triangle_popmap)
        tuple val(meta4), path(popmap)
    output:
        tuple val(meta), path("${meta.id}_triangle_mqc.html"), emit: plot_html
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    plot_triangle.py \\
        --result ${hindex} \\
        --result_fixed ${hindex_fixed} \\
        --popmap ${popmap} \\
        --triangle_map ${triangle_popmap} \\
        --template ${baseDir}/assets/multiqc_triangle.html \\
        --out "${meta.id}_triangle_mqc.html" \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
