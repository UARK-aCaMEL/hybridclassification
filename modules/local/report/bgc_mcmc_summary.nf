process BGC_MCMC_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta),        path(bgc_results)

    output:
        tuple val(meta), path("${meta.id}_bgc_summary_mqc.json"),      emit: table_json

    script:
    def args   = task.ext.args ?: ''

    """
    shopt -s nullglob

    gencline_matches=( "${bgc_results}"/*__gencline_hmc__summary.txt "${bgc_results}"/*__gencline__summary.txt )
    hindex_matches=( "${bgc_results}"/*__hi_hmc__summary.txt "${bgc_results}"/*__hindex_hmc__summary.txt "${bgc_results}"/*__hi__summary.txt )

    if (( \${#gencline_matches[@]} == 0 )); then
        echo "ERROR: Could not find gencline summary file in ${bgc_results}" >&2
        echo "Contents of ${bgc_results}:" >&2
        ls -lah "${bgc_results}" >&2 || true
        exit 1
    fi

    if (( \${#hindex_matches[@]} == 0 )); then
        echo "ERROR: Could not find hindex summary file in ${bgc_results}" >&2
        echo "Contents of ${bgc_results}:" >&2
        ls -lah "${bgc_results}" >&2 || true
        exit 1
    fi

    GENCLINE_SUMMARY="\${gencline_matches[0]}"
    HINDEX_SUMMARY="\${hindex_matches[0]}"

    echo "Using gencline summary: \${GENCLINE_SUMMARY}"
    echo "Using hindex summary: \${HINDEX_SUMMARY}"

    bgc_summary_table.py \\
        --gencline "\${GENCLINE_SUMMARY}" \\
        --hindex   "\${HINDEX_SUMMARY}" \\
        --template ${baseDir}/assets/multiqc_bgc_summary.html \\
        --out      "${meta.id}_bgc_summary_mqc.json" \\
        ${args}
    """
}
