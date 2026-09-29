//
// Subworkflow with functionality specific to the aCaMEL/hybridclassification pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { HTSLIB_BGZIPTABIX as BGZIP_INDEX_VCF } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             // string: Path to input VCF or VCF.gz file
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message
    popmap            // string: path to popmap file
    speciesmap
    site_coords
    geo_data_config
    geo_data_dir
    combinations

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Random seed for SNPio thinning, ADMIXTURE, hybrid simulation, NewHybrids and bgc
    //
    def random_seed = randomSeed()
    log.info(params.seed != null ? "Random seed: ${random_seed}" : "Random seed: ${random_seed} (not set; rerun with --seed ${random_seed} to reproduce this run)")

    //
    // Create channel from input file provided through params.input
    //

    Channel
        .fromPath(input)
        .map { file ->
            def meta = [id: file.simpleName]
            return [meta, file]
        }
        .branch {
            vcf: it[1].name.endsWith('.vcf')
            vcfgz: it[1].name.endsWith('.vcf.gz')
        }
        .set { ch_input }

    // Bgzip uncompressed VCFs (bgzipped inputs are linked as they are) and index them
    BGZIP_INDEX_VCF(
        ch_input.vcf.mix(ch_input.vcfgz).map { meta, file -> [ meta, file, [], [] ] },
        'compress',
        true,
        'vcf'
    )

    //
    // Create channel for the popmap
    //
    Channel
        .fromPath(popmap)
        .map { file ->
            def meta = [id: file.simpleName]
            return [meta, file]
        }
        .set{ ch_popmap }

    //
    // Channel for speciesmap
    //
    Channel
        .fromPath(speciesmap)
        .map { file ->
            def meta = [id: file.simpleName]
            return [meta, file]
        }
        .set{ ch_speciesmap }

    //
    // Channel for geo_data_config (optional)
    //
    if ( params.geo_data_config ) {
        Channel
            .fromPath( params.geo_data_config )
            .map { file ->
                def meta = [ id: file.simpleName ]
                return [ meta, file ]
            }
            .set { ch_geo_data_config }
    }
    else {
        Channel
            .empty()
            .set { ch_geo_data_config }
    }


    //
    // Channel for a *pre‑staged* geodata directory (optional)
    //
    if ( params.geo_data_dir ) {
        Channel
            .fromPath( params.geo_data_dir )      // accepts dir or wildcard
            .map { dir ->
                def meta = [ id: file(dir).getBaseName() ]
                return [ meta, dir ]
            }
            .set { ch_geo_data_dir }
    }
    else {
        Channel.empty().set { ch_geo_data_dir }
    }

    //
    // Channel for site_coords (optional)
    //
    if ( params.site_coords ) {
        Channel
            .fromPath( params.site_coords )
            .map { file ->
                def meta = [ id: file.simpleName ]
                return [ meta, file ]
            }
            .set { ch_site_coords }
    }
    else {
        Channel
            .empty()
            .set { ch_site_coords }
    }

    //
    // Channel for combinations to test
    //
    Channel
        .fromPath(combinations)
        .splitText()
        .map { line ->
            def (pop1, pop2) = line.tokenize()
            def meta = [
                id  : "${pop1}_${pop2}",
                pop1: pop1,
                pop2: pop2
            ]
            return [meta]
        }
        .set { ch_combinations }

    emit:
    vcf       = BGZIP_INDEX_VCF.out.output
    tbi       = BGZIP_INDEX_VCF.out.index
    popmap    = ch_popmap
    speciesmap = ch_speciesmap
    site_coords = ch_site_coords
    geo_data    = ch_geo_data_config
    geo_data_dir = ch_geo_data_dir
    combinations = ch_combinations
    seed      = channel.value(random_seed)
    versions  = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate input params
//
def validateInputParameters() {
    def errors = []

    def isIntegerLike = { v ->
        (v instanceof Integer) ||
        (v instanceof Long) ||
        (v instanceof Short) ||
        (v instanceof Byte)
    }

    def isNumericLike = { v ->
        (v instanceof Integer) ||
        (v instanceof Long) ||
        (v instanceof Short) ||
        (v instanceof Byte) ||
        (v instanceof Float) ||
        (v instanceof Double) ||
        (v instanceof BigDecimal)
    }

    def checkInteger = { name, value, min = null, max = null ->
        if (!isIntegerLike.call(value)) {
            errors << "Invalid value for --${name}: '${value}'. It must be an integer."
            return
        }
        if (min != null && value < min) {
            errors << "Invalid value for --${name}: '${value}'. It must be >= ${min}."
        }
        if (max != null && value > max) {
            errors << "Invalid value for --${name}: '${value}'. It must be <= ${max}."
        }
    }

    def checkNumeric = { name, value, min = null, max = null, minInclusive = true, maxInclusive = true ->
        if (!isNumericLike.call(value)) {
            errors << "Invalid value for --${name}: '${value}'. It must be numeric."
            return
        }

        if (min != null) {
            boolean badMin = minInclusive ? (value < min) : (value <= min)
            if (badMin) {
                String op = minInclusive ? ">=" : ">"
                errors << "Invalid value for --${name}: '${value}'. It must be ${op} ${min}."
            }
        }

        if (max != null) {
            boolean badMax = maxInclusive ? (value > max) : (value >= max)
            if (badMax) {
                String op = maxInclusive ? "<=" : "<"
                errors << "Invalid value for --${name}: '${value}'. It must be ${op} ${max}."
            }
        }
    }

    def checkBoolean = { name, value ->
        if (!(value instanceof Boolean)) {
            errors << "Invalid value for --${name}: '${value}'. It must be true or false."
        }
    }

    // Integer parameters
    checkInteger.call('maxk', params.maxk, 1)
    checkInteger.call('thin_dist', params.thin_dist, 0)
    checkInteger.call('sample_size', params.sample_size, 1)
    checkInteger.call('n_reps', params.n_reps, 1)
    checkInteger.call('panel_size', params.panel_size, 1)
    checkInteger.call('nh_burnin', params.nh_burnin, 0)
    checkInteger.call('nh_sweeps', params.nh_sweeps, 1)
    checkInteger.call('bgc_iters', params.bgc_iters, 1)
    checkInteger.call('bgc_thin', params.bgc_thin, 1)

    // Proportions / probabilities / bounded numeric parameters
    checkNumeric.call('ind_cov', params.ind_cov, 0, 1, true, true)
    checkNumeric.call('snp_cov', params.snp_cov, 0, 1, true, true)
    checkNumeric.call('pop_cov', params.pop_cov, 0, 1, true, true)
    checkNumeric.call('min_maf', params.min_maf, 0, 0.5, true, true)
    checkNumeric.call('ancestry_threshold', params.ancestry_threshold, 0, 1, true, true)
    checkNumeric.call('prob_threshold', params.prob_threshold, 0, 1, true, true)
    checkNumeric.call('af_dist_min', params.af_dist_min, 0, 1, true, true)
    checkNumeric.call('outlier_alpha', params.outlier_alpha, 0, 1, false, false)
    checkNumeric.call('bgc_burnin', params.bgc_burnin, 0, 1, false, true)

    // Boolean parameters
    checkBoolean.call('run_bgc', params.run_bgc)

    // Cross-parameter logic checks
    if (isIntegerLike.call(params.panel_size) && isIntegerLike.call(params.thin_dist)) {
        if (params.panel_size < 1) {
            errors << "Invalid value for --panel_size: '${params.panel_size}'. It must be >= 1."
        }
        if (params.thin_dist < 0) {
            errors << "Invalid value for --thin_dist: '${params.thin_dist}'. It must be >= 0."
        }
    }

    if (isIntegerLike.call(params.maxk) && params.maxk < 2) {
        log.warn "Parameter --maxk is '${params.maxk}'. Structure-like clustering usually expects K >= 2."
    }

    if (params.run_bgc instanceof Boolean && params.run_bgc) {
        if (isIntegerLike.call(params.bgc_iters) && isIntegerLike.call(params.bgc_thin)) {
            if (params.bgc_thin > params.bgc_iters) {
                errors << "Invalid combination: --bgc_thin (${params.bgc_thin}) cannot be greater than --bgc_iters (${params.bgc_iters})."
            }
        }
    }

    if (errors) {
        errors.each { log.error it }
        throw new IllegalArgumentException("Invalid input parameter(s). See error messages above.")
    }
}
//
// Random seed: --seed if given, otherwise derived from the session ID, which is
// new for every run but kept by -resume (so cached tasks stay valid)
//
def randomSeed() {
    return params.seed != null ? params.seed as long : Math.floorMod(workflow.sessionId.hashCode() as long, 2147483646L) + 1
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
