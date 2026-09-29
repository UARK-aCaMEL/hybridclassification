process EVALADMIX {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/mussmann/admixpipe:3.2.2'

    input:
    tuple val(meta), path(ped), path(map), path(pfiles), path(qfiles), path(qfiles_json), path(popmap), path(clumpak_output), path(major_clusters), path(cvruns_json), path(qfilepaths_json)

    output:
    tuple val(meta), path("*corres"), emit: corres
    tuple val(meta), path("[1-9]*.png"), emit: majorclust_png
    tuple val(meta), path("*MinClust*.png"), optional:true, emit: minorclust_png
    tuple val(meta), path("${meta.id}*.png"), emit: reps_png
    tuple val(meta), path("${meta.id}*.fam"), emit: fam
    path "versions.yml",  emit: versions

    script:
    def args   = task.ext.args ?: ''

    """
    # Dynamically add admixpipe paths if present in the container
    if [ -d /app ]; then
        export PATH="/app/bin:/app/scripts/python/clumpak:/app/scripts/python/admixturePipeline:\$PATH"
    fi

    runEvalAdmix.py \\
    -p ${meta.id} \\
    -k 1 \\
    -K ${params.maxk} \\
    -m ${popmap} \\
    -n ${task.cpus} \\
    ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AdmixPipe: 3.2.2
    END_VERSIONS
    """
}
