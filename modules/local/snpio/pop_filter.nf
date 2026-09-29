process SNPIO_POPFILTER {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/btmartin721/snpio:1.7.6'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(tbi)
    tuple val(meta3), path(popmap)

    output:
    tuple val(meta), path("${meta.id}.temp.filtered.vcf.gz"), emit: filtered_vcf
    tuple val(meta), path("${meta.id}.temp.filtered.vcf.gz.tbi"), emit: filtered_tbi
    tuple val(meta), path("*_output"), emit: snpio_output
    path "versions.yml",     emit: versions

    script:
    def args   = task.ext.args ?: ''

    """
    snpio_popfilter.py \\
        --vcf ${vcf} \\
        --popmap ${popmap} \\
        --pop_cov ${params.pop_cov} \\
        --prefix ${meta.id}".temp" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SNPio: \$(python -c "import snpio; print(snpio.__version__)")
    END_VERSIONS
    """
}
