process FIND_CANDIDATES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
        tuple val(meta),      path(k2_clumpp)
        tuple val(meta2),     path(inds)
        tuple val(meta3),     path(popmap)
        tuple val(meta4),     path(speciesmap)

    output:
        tuple val(meta), path("${meta.id}_popmap.tsv"), emit: popmap
        path "versions.yml"                    , emit: versions

    script:
    """
    gawk -v thresh=${params.ancestry_threshold} '\
    BEGIN { \
        OFS="\t"; idx=0; \
        while ((getline line < "${inds}") > 0) inds[++idx] = line; \
    } \
    { \
        sub(/.*: /, ""); \
        p1 = \$1; p2 = \$2; \
        maxval = (p1 > p2 ? p1 : p2); \
        if (maxval < thresh) assign = "ADMIX"; \
        else assign = (p1 > p2 ? "P0" : "P1"); \
        ind = inds[FNR]; \
        print ind, assign; \
    } \
    ' ${k2_clumpp} > ${meta.id}_popmap.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | grep -oP '(?<=GNU Awk ).*?(?=, )')
    END_VERSIONS
    """
}