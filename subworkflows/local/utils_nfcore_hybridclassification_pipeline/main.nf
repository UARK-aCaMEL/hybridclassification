//
// Subworkflow with functionality specific to the aCaMEL/hybridclassification pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIP               } from '../../../modules/nf-core/tabix/bgzip/main'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             // string: Path to input VCF or VCF.gz file
    popmap            // string: path to popmap file
    speciesmap
    site_coords
    geo_data_config
    geo_data_dir
    combinations

    main:

    ch_versions = Channel.empty()

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
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
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

    // Process VCF inputs
    TABIX_BGZIP ( ch_input.vcf )
    ch_tabix_vcf_input = ch_input.vcfgz
        | mix (TABIX_BGZIP.out.output )
    TABIX_TABIX( ch_tabix_vcf_input )

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

    // Collect versions
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
    vcf       = ch_tabix_vcf_input
    tbi       = TABIX_TABIX.out.tbi
    popmap    = ch_popmap
    speciesmap = ch_speciesmap
    site_coords = ch_site_coords
    geo_data    = ch_geo_data_config
    geo_data_dir = ch_geo_data_dir
    combinations = ch_combinations
    versions  = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
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
        if (!isIntegerLike(value)) {
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
        if (!isNumericLike(value)) {
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
    checkInteger('maxk', params.maxk, 1)
    checkInteger('thin_dist', params.thin_dist, 0)
    checkInteger('sample_size', params.sample_size, 1)
    checkInteger('n_reps', params.n_reps, 1)
    checkInteger('panel_size', params.panel_size, 1)
    checkInteger('nh_burnin', params.nh_burnin, 0)
    checkInteger('nh_sweeps', params.nh_sweeps, 1)
    checkInteger('bgc_iters', params.bgc_iters, 1)
    checkInteger('bgc_thin', params.bgc_thin, 1)

    // Proportions / probabilities / bounded numeric parameters
    checkNumeric('ind_cov', params.ind_cov, 0, 1, true, true)
    checkNumeric('snp_cov', params.snp_cov, 0, 1, true, true)
    checkNumeric('pop_cov', params.pop_cov, 0, 1, true, true)
    checkNumeric('min_maf', params.min_maf, 0, 0.5, true, true)
    checkNumeric('ancestry_threshold', params.ancestry_threshold, 0, 1, true, true)
    checkNumeric('prob_threshold', params.prob_threshold, 0, 1, true, true)
    checkNumeric('af_dist_min', params.af_dist_min, 0, 1, true, true)
    checkNumeric('outlier_alpha', params.outlier_alpha, 0, 1, false, false)
    checkNumeric('bgc_burnin', params.bgc_burnin, 0, 1, false, true)

    // Boolean parameters
    checkBoolean('run_bgc', params.run_bgc)

    // Cross-parameter logic checks
    if (isIntegerLike(params.panel_size) && isIntegerLike(params.thin_dist)) {
        if (params.panel_size < 1) {
            errors << "Invalid value for --panel_size: '${params.panel_size}'. It must be >= 1."
        }
        if (params.thin_dist < 0) {
            errors << "Invalid value for --thin_dist: '${params.thin_dist}'. It must be >= 0."
        }
    }

    if (isIntegerLike(params.maxk) && params.maxk < 2) {
        log.warn "Parameter --maxk is '${params.maxk}'. Structure-like clustering usually expects K >= 2."
    }

    if (params.run_bgc instanceof Boolean && params.run_bgc) {
        if (isIntegerLike(params.bgc_iters) && isIntegerLike(params.bgc_thin)) {
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
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
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
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
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
