process SNPIO_SELECT{
    tag "$pops.id"
    label 'process_medium'

    container 'docker.io/btmartin721/snpio:1.3.21'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(tbi)
    tuple val(meta3), path(popmap)
    val(pops)

    output:
    tuple val(pops), path("${pops.id}.subset.nremover.vcf.gz"), emit: filtered_vcf
    tuple val(pops), path("${pops.id}.subset.nremover.vcf.gz.tbi"), emit: filtered_tbi
    path "versions.yml",     emit: versions

    script:
    """
    snpio_subset.py \\
        --vcf ${vcf} \\
        --popmap ${popmap} \\
        --include ${pops.pop1} \\
        --include ${pops.pop2} \\
        --prefix ${pops.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SNPio: 1.3.21
    END_VERSIONS
    """
}
