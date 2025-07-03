process SIMULATE_HYBRIDS {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/tkchafin/hybriddetective:1.3'

    input:
    tuple val(meta), path(pure)

    output:
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    simulate_hybrids.R \\
        --panel ${pure} \\
        --pair "P0,P1" \\
        --n-sims ${params.n_sims} \\
        --n-reps ${params.n_reps} \\
        --prop-sample 2.0 \\
        --sample-size ${params.sample_size} \\
        --prefix ${meta.id} \\
        ${args}

    # grab hybriddetective version from R
    HYBRID_VER=\$(Rscript -e 'cat(as.character(utils::packageVersion("hybriddetective")))' )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hybriddetective: \$HYBRID_VER
    END_VERSIONS
    """
}
