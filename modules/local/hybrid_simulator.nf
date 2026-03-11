process SIMULATE_HYBRIDS {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/tkchafin/pysam:1.0'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(popmap)

    output:
    tuple val(meta), path("${meta.id}_simulation.vcf"),     emit: vcf
    tuple val(meta), path("${meta.id}_simulation.tsv"),     emit: popmap

    script:
    """
    hybrid_simulator.py \\
        --vcf ${vcf} \\
        --popmap ${popmap} \\
        --p0 "P0" \\
        --p1 "P1" \\
        --num_reps ${params.n_reps} \\
        --size_pure ${params.sample_size} \\
        --size_f1 ${params.sample_size} \\
        --size_f2 ${params.sample_size} \\
        --size_bc ${params.sample_size} \\
        --strategy "freq" \\
        --out_prefix ${meta.id}
    """
}