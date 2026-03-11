process VCF_TO_NEWHYBRIDS {
    tag   "$meta.id"
    label 'process_single'

    container 'docker.io/btmartin721/snpio:1.3.21'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(popmap)
    tuple val(meta3), path(sim_vcf)

    output:
    tuple val(meta), path("${meta.id}_NewHybrids.txt"),     emit: newhybrids
    tuple val(meta), path("${meta.id}_NH_index_map.tsv"),   emit: index_map

    script:
    """
    # build --simulation flag only if we got a sim_vcf input
    SIM_ARG=""
    if [[ -n "${sim_vcf}" ]]; then
        SIM_ARG="--simulation ${sim_vcf}"
    fi

    vcf_to_genepop.py \
        --vcf        ${vcf} \
        --popmap     ${popmap} \
        --prefix     ${meta.id} \
        --p0-label   P0 \
        --p1-label   P1 \
        \$SIM_ARG
    """
}
