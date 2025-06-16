process INPUT_TABLE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
        tuple val(meta), path(k2_clumpp)
        tuple val(meta2), path(inds)
        tuple val(meta3), path(popmap)
        tuple val(meta4), path(speciesmap)

    output:
        tuple val(meta), path("${meta.id}_assignments.tsv")       , emit: table
        path "versions.yml"                                       , emit: versions

    script:
    """


    # Log versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | grep -oP '(?<=GNU Awk ).*?(?=, )')
    END_VERSIONS
    """
}
