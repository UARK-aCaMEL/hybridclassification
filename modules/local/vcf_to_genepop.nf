process VCF_TO_GENEPOP {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/btmartin721/snpio:1.3.21'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(popmap)

    output:
    tuple val(meta), path("${meta.id}_pure.txt"), emit: pure
    tuple val(meta), path("${meta.id}_admix.txt"), emit: admix

    script:

    """
    vcf_to_genepop.py \\
        ${vcf} \\
        ${popmap} \\
        ${meta.id}
    """
}
