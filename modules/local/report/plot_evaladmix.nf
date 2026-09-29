process PLOT_EVALADMIX {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(qfilepaths), path(fam), path(corres)
        path(template)

    output:
        tuple val(meta), path("${meta.id}_evaladmix_k2_mqc.html"), emit: plot_html
        path("versions.yml")   , emit: versions

    script:
    """
    # Candidate hybrids come from K = 2, so plot the fit of that model
    plot_evaladmix.py \\
            --prefix "${meta.id}" \\
            --workdir . \\
            --qfilePaths ${qfilepaths} \\
            --best_k 2 \\
            --bestk_html "${meta.id}_evaladmix_k2_mqc.html" \\
            --template_bestk ${template}

    plotly_version=\$(python3 -c 'import plotly; print(plotly.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
