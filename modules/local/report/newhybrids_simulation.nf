process NH_PLOT_SIMULATION {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(sim_result)
        tuple val(meta2), path(sim_map)
    output:
        tuple val(meta), path("${meta.id}_nh_sim_mqc.html"), emit: plot_html
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    plot_nh_simulation.py \\
        --result ${sim_result} \\
        --map ${sim_map} \\
        --threshold ${params.prob_threshold} \\
        --template ${baseDir}/assets/multiqc_nh_sim.html \\
        --out "${meta.id}_nh_sim_mqc.html" \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
