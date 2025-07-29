process NH_SUMMARY_TABLE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta),        path(nh_results)
        tuple val(meta2),       path(nh_map)
        tuple val(meta3),       path(popmap)
        tuple val(meta4),       path(speciesmap)
        tuple val(meta5),       path(masked_samples)

    output:
        tuple val(meta), path("${meta.id}_nh_summary_mqc.json"),         emit: table_json
        tuple val(meta), path("${meta.id}_nh_hybrids.txt"),             emit: hybrid_list
        tuple val(meta), path("${meta.id}_nh_hybrids_masked.txt"),      emit: hybrid_list_masked,   optional: true
        tuple val(meta), path("${meta.id}_nh_summary_mask_mqc.json"),   emit: table_json_masked,    optional: true
        path("versions.yml"),                                          emit: versions

    script:
    def args   = task.ext.args ?: ''
    // build the mask arguments only if we actually got a masked_samples file
    def maskArgs = masked_samples ?
        " --mask ${masked_samples} --masked_template ${baseDir}/assets/multiqc_nh_summary_masked.html --out_mask ${meta.id}_nh_summary_mask_mqc.json --list_masked ${meta.id}_nh_hybrids_masked.txt "
      : ""

    """
    nh_summary_table.py \\
        --result       ${nh_results} \\
        --result_map   ${nh_map} \\
        --popmap       ${popmap} \\
        --speciesmap   ${speciesmap} \\
        --threshold    ${params.prob_threshold} \\
        --template     ${baseDir}/assets/multiqc_nh_summary.html \\
        --out          \"${meta.id}_nh_summary_mqc.json\" \\
        --list         \"${meta.id}_nh_hybrids.txt\" \\
        ${maskArgs} \\
        ${args}

    plotly_version=\$(python3 - <<EOF
    import plotly
    print(plotly.__version__)
    EOF
    )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \${plotly_version}
    END_VERSIONS
    """
}
