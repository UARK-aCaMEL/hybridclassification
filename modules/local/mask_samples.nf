process MASK_SAMPLES {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/btmartin721/snpio:1.3.21'

    input:
    tuple val(meta), path(hindex)
    tuple val(meta2), path(hindex_popmap)
    tuple val(meta3), path(pofz)
    tuple val(meta3), path(index_map)

    output:
    tuple val(meta), path("${meta.id}_masked_samples.txt"),     emit: mask
    tuple val(meta), path("${meta.id}_mask_table.tsv"),     emit: mask_table

    script:
    def args = task.ext.args ?: ''
    """
    mask_samples.py \\
        --hindex ${hindex} \\
        --nh_results ${pofz} \\
        --nh_index ${index_map} \\
        --hindex_popmap ${hindex_popmap} \\
        --out_prefix ${meta.id} \\
        --prob_threshold ${params.prob_threshold} \\
        --alpha ${params.outlier_alpha} \\
        ${args}
    """
}