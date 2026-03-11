process NH_PLOT_TRACE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(nh_trace)
    output:
        tuple val(meta), path("${meta.id}_nh_trace_mqc.html"), emit: plot_html
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    plot_nh_trace.py \\
        --trace ${nh_trace} \\
        --template ${baseDir}/assets/multiqc_nh_trace.html \\
        --burnin ${params.nh_burnin} \\
        --out "${meta.id}_nh_trace_mqc.html" \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
