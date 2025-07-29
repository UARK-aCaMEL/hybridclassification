process TRIANGLE_METRICS {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/tkchafin/pysam:1.0'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(popmap)
    tuple val(meta3), path(sim_vcf)

    output:
    tuple val(meta), path("${meta.id}_classification_popmap.tsv"),     emit: popmap
    tuple val(meta), path("${meta.id}_hindex.tsv"),     emit: hindex
    tuple val(meta), path("${meta.id}_hindex_fixed.tsv"),     emit: hindex_fixed

    script:
    def args = task.ext.args ?: ''
    """
    # build --simulation flag only if we received a sim_vcf input
    SIM_ARG=""
    if [[ -n "${sim_vcf}" ]]; then
        SIM_ARG="--simulation ${sim_vcf}"
    fi
    triangle_metrics.py \\
        --vcf ${vcf} \\
        --popmap ${popmap} \\
        --p0 "P0" \\
        --p1 "P1" \\
        --out_prefix ${meta.id} \\
        --af_diff_min ${params.af_dist_min} \\
        \$SIM_ARG \\
        ${args}
    """
}