process BGC_PLOT_HINDEX {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/tkchafin/plotly:1.1"

    input:
        tuple val(meta), path(bgc_results)

    output:
        tuple val(meta), path("${meta.id}_hindex_plot_mqc.html"), emit: plot_html

    script:
    def args = task.ext.args ?: ''

    """
    shopt -s nullglob

    order_matches=( "${bgc_results}"/*__input_order__GlikADMIX_samples.txt )
    hi_matches=( "${bgc_results}"/*__h_out__hi.tsv )

    if (( \${#order_matches[@]} == 0 )); then
        echo "ERROR: Could not find hybrid index sample order file in ${bgc_results}" >&2
        echo "Contents of ${bgc_results}:" >&2
        ls -lah "${bgc_results}" >&2 || true
        exit 1
    fi

    if (( \${#hi_matches[@]} == 0 )); then
        echo "ERROR: Could not find hybrid index results file in ${bgc_results}" >&2
        echo "Contents of ${bgc_results}:" >&2
        ls -lah "${bgc_results}" >&2 || true
        exit 1
    fi

    ORDER_FILE="\${order_matches[0]}"
    HI_FILE="\${hi_matches[0]}"

    echo "Using input order file: \${ORDER_FILE}"
    echo "Using hybrid index results: \${HI_FILE}"

    bgc_plot_hindex.py \\
        --order    "\${ORDER_FILE}" \\
        --hindex   "\${HI_FILE}" \\
        --template ${baseDir}/assets/multiqc_bgc_hindex.html \\
        --out      "${meta.id}_hindex_plot_mqc.html" \\
        ${args}
    """
}
