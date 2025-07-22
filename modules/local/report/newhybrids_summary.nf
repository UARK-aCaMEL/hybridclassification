process NH_SUMMARY_TABLE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(nh_results)
        tuple val(meta2), path(nh_map)
        tuple val(meta3), path(popmap)
        tuple val(meta4), path(speciesmap)
    output:
        tuple val(meta), path("${meta.id}_nh_summary_mqc.json"), emit: table_json
        tuple val(meta), path("${meta.id}_nh_hybrids.txt"), emit: hybrid_list
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    nh_summary_table.py \\
        --result ${nh_results} \\
        --result_map ${nh_map} \\
        --popmap ${popmap} \\
        --speciesmap ${speciesmap} \\
        --threshold ${params.prob_threshold} \\
        --template ${baseDir}/assets/multiqc_nh_summary.html \\
        --out "${meta.id}_nh_summary_mqc.json" \\
        --list "${meta.id}_nh_hybrids.txt" \\
        ${args}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
