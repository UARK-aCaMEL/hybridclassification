process BESTK {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
        tuple val(meta), path(cv_file)
        tuple val(meta2), path(best_results)

    output:
        tuple val(meta), path('bestK.txt')                        , emit: bestK_file
        tuple val(meta), path('best_clumpp_indfile.out')          , emit: bestK_clumpp
        tuple val(meta), path('k2_clumpp_indfile.out')          , emit: k2_clumpp
        path "versions.yml"                                       , emit: versions

    script:
    """
    # Determine best K (excluding K=1)
    awk 'NR>1 && \$1!=1 {
        if (min=="" || \$2<min) { min=\$2; k=\$1 }
    }
    END { print k }' ${cv_file} > bestK.txt

    # Read best K
    K=\$(cat bestK.txt)

    # Locate the matching ClumppIndFile.output.\$K
    MATCHED_FILE=\$(find "${best_results}/" -type f -name "ClumppIndFile.output.\$K" -print)

    if [ ! -f "\$MATCHED_FILE" ]; then
        echo "❌ Could not find ClumppIndFile.output.\$K in ${best_results}" >&2
        exit 1
    fi

    # Symlink it to a standard name for downstream use
    MATCHED_FILE=\$(realpath "\$MATCHED_FILE")
    ln -s "\$MATCHED_FILE" best_clumpp_indfile.out

    # Also locate the file for K=2
    K2_MATCH=\$(find "${best_results}/" -type f -name "ClumppIndFile.output.2" -print)
    if [ ! -f "\$K2_MATCH" ]; then
        echo "❌ Could not find ClumppIndFile.output.2 in ${best_results}" >&2
        exit 1
    fi
    ln -s "\$(realpath "\$K2_MATCH")" k2_clumpp_indfile.out

    # Log versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | grep -oP '(?<=GNU Awk ).*?(?=, )')
    END_VERSIONS
    """
}
