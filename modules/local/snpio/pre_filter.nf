process SNPIO_FILTER {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/btmartin721/snpio:1.3.6'

    input:
    tuple val(meta), path(vcf), path(popmap)

    output:
    tuple val(meta), path("${meta.id}.filter.nremover.vcf.gz"), emit: filtered_vcf
    tuple val(meta), path("${meta.id}.filter.nremover.vcf.gz.tbi"), emit: filtered_tbi
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
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SNPio: 1.3.6
    END_VERSIONS
    """
}
