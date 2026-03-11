process NH_PLOT_SPATIAL {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/geopandas:1.0"

    input:
        tuple val(meta), path(nh_results)
        tuple val(meta2), path(nh_map)
        tuple val(meta3), path(popmap)
        tuple val(meta4), path(site_coords)
        tuple val(meta5), path(geo_data)
        tuple val(meta6), path(mask)
    output:
        tuple val(meta), path("${meta.id}_nh_spatial_mqc.html"), emit: plot_html
        tuple val(meta), path("${meta.id}_points.tsv"), emit: table
        path("versions.yml")   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def geo_data_arg = geo_data ? "--geo_data_json ${geo_data}/config.json" : ''
    """
    plot_nh_spatial.py \\
        --result ${nh_results} \\
        --result_map ${nh_map} \\
        --popmap ${popmap} \\
        --template ${baseDir}/assets/multiqc_nh_spatial.html \\
        --out "${meta.id}_nh_spatial_mqc.html" \\
        --site_coords ${site_coords} \\
        --threshold ${params.prob_threshold} \\
        --table_out "${meta.id}_points.tsv" \\
        --mask ${mask} \\
        ${geo_data_arg} \\
        ${args}

    geopandas_version=\$(python3 -c 'import geopandas; print(geopandas.__version__)')
    folium_version=\$(python3 -c 'import folium; print(folium.__version__)')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        geopandas: \${geopandas_version}
        folium: \${folium_version}
    END_VERSIONS
    """
}
