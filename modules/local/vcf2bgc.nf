process VCF2BGC {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/tkchafin/pysam:1.0'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(popmap)

    output:
    tuple val(meta), path("${meta.id}_GlikP0.txt"),            emit: parent1
    tuple val(meta), path("${meta.id}_GlikP1.txt"),            emit: parent2
    tuple val(meta), path("${meta.id}_GlikADMIX.txt"),         emit: admix

    tuple val(meta), path("${meta.id}_GlikP0_samples.txt"),    emit: parent1_samples
    tuple val(meta), path("${meta.id}_GlikP1_samples.txt"),    emit: parent2_samples
    tuple val(meta), path("${meta.id}_GlikADMIX_samples.txt"), emit: admix_samples
    tuple val(meta), path("${meta.id}_locus_order.txt"),       emit: locus_order

    script:
    def args = task.ext.args ?: ''
    """
    vcf2bgc.py \\
        --vcf ${vcf} \\
        --popmap ${popmap} \\
        --p0 "P0" \\
        --p1 "P1" \\
        --out_prefix ${meta.id} \\
        --af_diff_min ${params.af_dist_min} \\
        ${args}
    """
}
