process NH_PLOT_CLASSIFICATIONS {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(nh_results)
        tuple val(meta2), path(nh_map)
        tuple val(meta3), path(popmap)
        tuple val(meta4), path(speciesmap)
    output:
        tuple val(meta), path("${meta.id}_nh_classifications_mqc.html"), emit: plot_html
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    plot_nh_classifications.py \\
        --result ${nh_results} \\
        --result_map ${nh_map} \\
        --popmap ${popmap} \\
        --speciesmap ${speciesmap} \\
        --template ${baseDir}/assets/multiqc_nh_classifications.html \\
        --out "${meta.id}_nh_classifications_mqc.html" \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
