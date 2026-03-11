process RUN_NEWHYBRIDS {
    tag   "$meta.id"
    label 'process_medium'
    container 'docker.io/tkchafin/hybriddetective:1.4'

    input:
    tuple val(meta), path(nh_file)

    output:
    path "versions.yml",                     emit: versions
    tuple val(meta), path("${meta.id}_nh_results.txt"),        emit: output
    tuple val(meta), path("EchoedGtypData.txt"),               emit: echoedData
    tuple val(meta), path("aa-EchoedGtypFreqCats.txt"),        emit: freqCats
    tuple val(meta), path("aa-LociAndAlleles.txt"),            emit: lociAlleles
    tuple val(meta), path("aa-Pi.aves"),                       emit: piAverages
    tuple val(meta), path("aa-Pi.hist"),                       emit: piHist
    tuple val(meta), path("aa-PofZ.txt"),                      emit: pofz
    tuple val(meta), path("aa-ScaledLikelihood.txt"),          emit: scaledLikelihood
    tuple val(meta), path("aa-Theta.hist"),                    emit: thetaHist
    tuple val(meta), path("aa-ThetaAverages.txt"),             emit: thetaAverages
    tuple val(meta), path("pi_trace.tsv"),                     emit: pi_trace

    script:
    def maxInt = (2**31) - 1
    def randOffset = new Random().nextInt(maxInt)
    def seed1 = Math.abs((task.index.toInteger() + randOffset) % maxInt)
    if( seed1 == 0 ) seed1 = 1
    def seed2 = (seed1 + 1) % maxInt
    if( seed2 == 0 ) seed2 = 1

    def args = task.ext.args ?: ''
    """
    # Because newhybs for some reason keeps giving a non-zero exit code
    set +e

    # run NewHybrids
    newhybrids \\
        --data-file    ${nh_file} \\
        --burn-in      ${params.nh_burnin} \\
        --num-sweeps   ${params.nh_sweeps} \\
        --seeds        ${seed1} ${seed2} \\
        --no-gui \\
        --print-traces Pi 1 \\
        --pi-prior "uniform" \\
        ${args} > ${meta.id}_nh_results.txt  2>&1

    # Grep out the Pi trace
    grep "PI_TRACE" ${meta.id}_nh_results.txt > pi_trace.tsv

    # record version - quote the output to handle special characters
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        newhybrids: "\$(newhybrids --version 2>&1 | head -n1)"
    END_VERSIONS

    # require that NewHybrids completed properly
    grep -q 'Output is in the following files:' ${meta.id}_nh_results.txt
    exit \$?
    """
}
