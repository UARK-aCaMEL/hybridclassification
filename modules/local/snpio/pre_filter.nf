process SNPIO_FILTER {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/btmartin721/snpio:1.7.6'

    input:
    tuple val(meta), path(vcf), path(popmap)
    val(seed)

    output:
    tuple val(meta), path("${meta.id}.filter.vcf.gz"), emit: filtered_vcf
    tuple val(meta), path("${meta.id}.filter.vcf.gz.tbi"), emit: filtered_tbi
    tuple val(meta), path("*_output"), emit: snpio_output
    path "versions.yml",     emit: versions

    script:
    def args   = task.ext.args ?: ''

    """
    snpio_filter.py \\
        --vcf ${vcf} \\
        --popmap ${popmap} \\
        --ind_cov ${params.ind_cov} \\
        --flank_dist ${params.thin_dist} \\
        --min_maf ${params.min_maf} \\
        --snp_cov ${params.snp_cov} \\
        --seed ${seed} \\
        --prefix ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SNPio: \$(python -c "import snpio; print(snpio.__version__)")
    END_VERSIONS
    """
}
