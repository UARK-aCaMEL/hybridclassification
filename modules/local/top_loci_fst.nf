process TOP_LOCI_FST {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.16--he513fc3_4' :
        'biocontainers/vcftools:0.1.16--he513fc3_4' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(tbi)
    tuple val(meta3), path(popmap)

    output:
    tuple val(meta), path("${meta.id}_top${topN}_fst.vcf.gz"), emit: top_vcf
    tuple val(meta), path("${meta.id}_fst_per_locus.tsv"),     emit: fst_table

    script:
    """
    # split popmap
    awk '\$2=="P0"{print \$1}' ${popmap} > P0.ids
    awk '\$2=="P1"{print \$1}' ${popmap} > P1.ids

    # compute Fst
    vcftools \\
    --gzvcf ${vcf} \\
    --weir-fst-pop P0.ids \\
    --weir-fst-pop P1.ids \\
    --out ${meta.id}_fst

    # select top N
    ( head -n1 ${meta.id}_fst.weir.fst && \\
    tail -n+2 ${meta.id}_fst.weir.fst \\
        | sort -k4,4nr \\
        | head -n ${params.panel_size} ) > ${meta.id}_fst_per_locus.tsv

    # subset VCF
    cut -f1,2 ${meta.id}_fst_per_locus.tsv > positions.txt

    vcftools \\
    --gzvcf ${vcf} \\
    --positions positions.txt \\
    --recode \\
    --stdout \\
    | bgzip -c > ${meta.id}_top${topN}_fst.vcf.gz

    tabix -p vcf ${meta.id}_top${topN}_fst.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(echo \$(vcftools --version 2>&1) | sed 's/^.*VCFtools (//;s/).*//')
    END_VERSIONS
    """
}
