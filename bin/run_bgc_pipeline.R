#!/usr/bin/env Rscript

# ------------------------------------------------------------
# BGC / genomic cline pipeline runner (3-input version)
# Inputs (dget):
#   - GlikP0    : parental population 0
#   - GlikP1    : parental population 1
#   - GlikADMIX : admixed/hybrid individuals (ONLY admix set)
#
# Expected genotype-likelihood format (glik):
#   Each object is typically a LIST of 4 nucleotide arrays (A,C,G,T),
#   where each array is an Ind x Loc matrix (Loc = columns).
#
# Runs:
#   est_p -> est_hi -> est_genocl -> sum2zero -> gencline_plot -> est_Q -> tri_plot
#
# Outputs:
#   outdir/
#     plots/  : PDF plots
#     text/   : text exports of sz_out (+ key outputs)
#     rdata/  : sz_out saved as .rds and .RData (+ all_outputs.rds)
# ------------------------------------------------------------

# ---- helpers ----
stop2 <- function(...) {
    message(...)
    quit(status = 1)
}

get_argval <- function(args, key, default = NULL) {
    # Supports: --key value  OR  --key=value
    hit <- grep(paste0("^", key, "($|=)"), args)
    if (length(hit) == 0) {
        return(default)
    }
    token <- args[hit[1]]
    if (grepl("=", token)) {
        sub(paste0("^", key, "="), "", token)
    } else {
        if (hit[1] == length(args)) stop2("Missing value after ", key)
        args[hit[1] + 1]
    }
}

ensure_dir <- function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    normalizePath(d, winslash = "/", mustWork = FALSE)
}

safe_dget <- function(path, label) {
    if (is.null(path) || is.na(path) || !nzchar(path)) stop2("Missing required path for ", label)
    if (!file.exists(path)) stop2("File not found for ", label, ": ", path)
    dget(path)
}

save_plot_pdf <- function(filename, expr, width = 7, height = 5) {
    pdf(filename, width = width, height = height, onefile = TRUE)
    on.exit(dev.off(), add = TRUE)
    force(expr)
}

write_any <- function(x, file_base, outdir) {
    # Writes object to sensible text formats:
    # - list (non-data.frame) -> recurse
    # - matrix/data.frame -> TSV (with rownames)
    # - atomic/scalar/other -> TXT (one per line or str() fallback)

    sanitize <- function(s) {
        s <- gsub("[^A-Za-z0-9._-]+", "_", s)
        s <- gsub("_+", "_", s)
        s
    }

    if (is.list(x) && !is.data.frame(x)) {
        nms <- names(x)
        if (is.null(nms)) nms <- paste0("elt", seq_along(x))
        for (i in seq_along(x)) {
            nm <- sanitize(nms[i])
            write_any(x[[i]], paste0(file_base, "__", nm), outdir)
        }
        return(invisible(NULL))
    }

    if (is.data.frame(x) || is.matrix(x)) {
        fn <- file.path(outdir, paste0(file_base, ".tsv"))
        write.table(
            x,
            file = fn, sep = "\t", quote = FALSE, col.names = NA,
            row.names = TRUE
        )
        return(invisible(NULL))
    }

    fn <- file.path(outdir, paste0(file_base, ".txt"))
    con <- file(fn, open = "wt")
    on.exit(close(con), add = TRUE)

    if (length(x) == 0) {
        writeLines("<empty>", con)
    } else if (is.atomic(x)) {
        writeLines(as.character(x), con)
    } else {
        writeLines(capture.output(str(x)), con)
    }
    invisible(NULL)
}

# ---- format checks ----

extract_glik_prefix <- function(path) {
    # Return the part before "_Glik" from the basename.
    # If not found, return NA.
    b <- basename(path)
    m <- regexpr("_Glik", b, fixed = TRUE)
    if (m[1] == -1) {
        return(NA_character_)
    }
    substr(b, 1, m[1] - 1)
}

infer_loc_n_from_matrix <- function(mat) {
    if (!is.matrix(mat) && !is.data.frame(mat)) {
        return(NA_integer_)
    }
    as.integer(ncol(mat))
}

get_glik_component_mats <- function(G) {
    # Returns list of matrices (A/C/G/T etc.) for validation.
    #
    # Typical: list(A=mat, C=mat, G=mat, T=mat) OR unnamed list of matrices.
    # We'll keep only matrix/data.frame elements.
    if (!is.list(G) || is.data.frame(G) || is.matrix(G)) {
        return(list())
    }
    mats <- G[vapply(G, function(x) is.matrix(x) || is.data.frame(x), logical(1))]
    mats
}

validate_glik_loci <- function(G, label) {
    mats <- get_glik_component_mats(G)

    if (length(mats) == 0) {
        stop2(
            label, " does not look like a glik list of nucleotide matrices.\n",
            "Expected a list where each element is an Ind x Loc matrix (Loc = columns)."
        )
    }

    loc_counts <- vapply(mats, infer_loc_n_from_matrix, integer(1))
    if (any(is.na(loc_counts))) {
        stop2(label, " contains non-matrix elements where matrices were expected.")
    }

    uniq <- unique(loc_counts)
    if (length(uniq) != 1) {
        nm <- names(mats)
        if (is.null(nm)) nm <- paste0("elt", seq_along(mats))
        msg <- paste0(nm, "=", loc_counts, collapse = ", ")
        stop2(
            label, " nucleotide matrices do not share the same number of loci (columns).\n",
            "Found: ", msg
        )
    }

    list(n_loci = uniq[1], component_names = names(mats))
}


write_stan_draws_wide <- function(fit, file_base, outdir) {
    if (is.null(fit)) {
        return(invisible(NULL))
    }

    df <- as.data.frame(fit)

    fn <- file.path(outdir, paste0(file_base, "__draws_wide.tsv"))
    write.table(
        df,
        file = fn,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    invisible(df)
}

write_stan_draws_long <- function(fit, file_base, outdir) {
    if (is.null(fit)) {
        return(invisible(NULL))
    }

    wide <- as.data.frame(fit)

    n_draws <- nrow(wide)
    n_chains <- fit@sim$chains
    n_iter_per_chain <- n_draws / n_chains

    if (n_iter_per_chain != floor(n_iter_per_chain)) {
        stop2("Could not evenly divide posterior draws by chains for ", file_base)
    }

    n_iter_per_chain <- as.integer(n_iter_per_chain)

    wide$chain <- rep(seq_len(n_chains), each = n_iter_per_chain)
    wide$iteration <- rep(seq_len(n_iter_per_chain), times = n_chains)

    param_cols <- setdiff(names(wide), c("chain", "iteration"))

    long_list <- lapply(param_cols, function(nm) {
        data.frame(
            parameter = nm,
            chain = wide$chain,
            iteration = wide$iteration,
            value = wide[[nm]],
            stringsAsFactors = FALSE
        )
    })

    long_df <- do.call(rbind, long_list)

    fn <- file.path(outdir, paste0(file_base, "__draws_long.tsv"))
    write.table(
        long_df,
        file = fn,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    invisible(long_df)
}

write_stan_summary <- function(fit, file_base, outdir) {
    if (is.null(fit)) {
        return(invisible(NULL))
    }

    txt <- capture.output(print(fit))

    fn <- file.path(outdir, paste0(file_base, "__summary.txt"))
    writeLines(txt, con = fn)

    invisible(txt)
}

write_stan_sampler_diagnostics <- function(fit, file_base, outdir, inc_warmup = FALSE) {
    if (is.null(fit)) {
        return(invisible(NULL))
    }

    sp <- rstan::get_sampler_params(fit, inc_warmup = inc_warmup)

    out_list <- lapply(seq_along(sp), function(ch) {
        df <- as.data.frame(sp[[ch]])
        df$chain <- ch
        df$iteration <- seq_len(nrow(df))
        df
    })

    diag_df <- do.call(rbind, out_list)

    suffix <- if (inc_warmup) {
        "__sampler_diagnostics_with_warmup.tsv"
    } else {
        "__sampler_diagnostics.tsv"
    }

    fn <- file.path(outdir, paste0(file_base, suffix))
    write.table(
        diag_df,
        file = fn,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    invisible(diag_df)
}

# ---- parse CLI args ----
args <- commandArgs(trailingOnly = TRUE)

if (any(args %in% c("-h", "--help"))) {
    cat(
        "Usage:
  run_bgc_pipeline.R --p0 PATH --p1 PATH --admix PATH
                    [--outdir DIR] [--prefix NAME]
                    [--n_iters 4000] [--warmup 2000] [--n_thin 1] [--ci 0.90]

Inputs are read with dget().

Required:
  --p0       Path to GlikP0 dget() file
  --p1       Path to GlikP1 dget() file
  --admix    Path to GlikADMIX dget() file (admixed individuals; used as Gx)

Optional:
  --outdir   Output directory (default: ./bgc_out)
  --prefix   Output prefix (default: derived from admix basename)
  --n_iters  Iterations for est_hi / est_genocl where supported (default: 4000)
  --warmup   Warmup proportion for est_hi / est_genocl where supported (default: 0.5)
  --n_thin   Thinning interval for est_hi / est_genocl where supported (default: 1)
  --ci       Credible interval level for sum2zero (default: 0.90)

Added validations:
  1) P0/P1/ADMIX glik objects must each have consistent loci counts across their nucleotide matrices,
     and all three must match each other.
  2) Input filenames must share the same prefix before \"_Glik\".
"
    )
    quit(status = 0)
}

p0_path <- get_argval(args, "--p0")
p1_path <- get_argval(args, "--p1")
admix_path <- get_argval(args, "--admix")

# ---- filename prefix check ----
p0_pref <- extract_glik_prefix(p0_path)
p1_pref <- extract_glik_prefix(p1_path)
ax_pref <- extract_glik_prefix(admix_path)

if (any(is.na(c(p0_pref, p1_pref, ax_pref)))) {
    stop2(
        "All input files must include \"_Glik\" in the basename so a prefix can be inferred.\n",
        "Basenames:\n",
        "  p0:    ", basename(p0_path), "\n",
        "  p1:    ", basename(p1_path), "\n",
        "  admix: ", basename(admix_path)
    )
}

if (length(unique(c(p0_pref, p1_pref, ax_pref))) != 1) {
    stop2(
        "Input files do not share the same prefix before \"_Glik\".\n",
        "Found prefixes:\n",
        "  p0:    ", p0_pref, "\n",
        "  p1:    ", p1_pref, "\n",
        "  admix: ", ax_pref, "\n",
        "Example expected: CAMANO_CAMOLI_GlikP0.txt / CAMANO_CAMOLI_GlikP1.txt / CAMANO_CAMOLI_GlikADMIX.txt"
    )
}

outdir <- ensure_dir(get_argval(args, "--outdir", default = "bgc_out"))

prefix <- get_argval(args, "--prefix", default = NA_character_)
if (is.na(prefix) || !nzchar(prefix)) {
    # Prefer inferred prefix from filenames when available
    prefix <- ax_pref
}

n_iters <- as.integer(get_argval(args, "--n_iters", default = "4000"))
warmup <- as.numeric(get_argval(args, "--warmup", default = "0.5"))
n_thin <- as.integer(get_argval(args, "--n_thin", default = "1"))
ci <- as.numeric(get_argval(args, "--ci", default = "0.90"))

if (is.na(n_iters) || n_iters <= 0) stop2("--n_iters must be a positive integer.")
if (is.na(warmup) || warmup < 0 || warmup >= 1) stop2("--warmup must be a numeric proportion in [0, 1).")
if (is.na(n_thin) || n_thin <= 0) stop2("--n_thin must be a positive integer.")
if (is.na(ci) || ci <= 0 || ci >= 1) stop2("--ci must be in (0, 1).")

plots_dir <- ensure_dir(file.path(outdir, "plots"))
text_dir <- ensure_dir(file.path(outdir, "text"))
rds_dir <- ensure_dir(file.path(outdir, "rdata"))

# ---- load package providing required functions ----
pkg_candidates <- c("bgc", "bgcTools", "bgchm", "bgc_hm", "bgc.hm", "bgc-hm")

pkg_loaded <- FALSE
for (p in pkg_candidates) {
    if (requireNamespace(p, quietly = TRUE)) {
        suppressPackageStartupMessages(library(p, character.only = TRUE))
        pkg_loaded <- TRUE
        message("Loaded package: ", p)
        break
    }
}
if (!pkg_loaded) {
    stop2(
        "Could not find a package providing: est_p, est_hi, est_genocl, sum2zero, gencline_plot, est_Q, tri_plot.\n",
        "Install/load the correct BGC package and/or edit pkg_candidates in this script."
    )
}

# ---- read inputs ----
GlikP0 <- safe_dget(p0_path, "GlikP0 (--p0)")
GlikP1 <- safe_dget(p1_path, "GlikP1 (--p1)")
GlikADMIX <- safe_dget(admix_path, "GlikADMIX (--admix)")

# ---- loci count checks ----
v0 <- validate_glik_loci(GlikP0, "GlikP0")
v1 <- validate_glik_loci(GlikP1, "GlikP1")
vx <- validate_glik_loci(GlikADMIX, "GlikADMIX")

if (length(unique(c(v0$n_loci, v1$n_loci, vx$n_loci))) != 1) {
    stop2(
        "Input datasets do not have the same number of loci (columns).\n",
        "  GlikP0:    ", v0$n_loci, "\n",
        "  GlikP1:    ", v1$n_loci, "\n",
        "  GlikADMIX: ", vx$n_loci
    )
}
message("Loci count check OK: ", v0$n_loci, " loci in P0/P1/ADMIX.")

# In downstream calls, the admixed set is used as Gx; keep your naming convention:
GlikHybrids <- GlikADMIX

# ---- analyses ----
message("Estimating parental allele frequencies (est_p) ...")
p_out <- est_p(G0 = GlikP0, G1 = GlikP1, model = "glik", ploidy = "diploid", HMC = FALSE)

message("Estimating hybrid indexes (est_hi) ...")
h_out <- est_hi(
    Gx = GlikHybrids,
    p0 = p_out$p0[, 1],
    p1 = p_out$p1[, 1],
    model = "glik",
    ploidy = "diploid",
    n_iters = n_iters,
    p_warmup = warmup,
    n_thin = n_thin
)

# ---- plot: hybrid index with intervals ----
hi_pdf <- file.path(plots_dir, paste0(prefix, "_hybrid_index.pdf"))
save_plot_pdf(hi_pdf, {
    ord <- order(h_out$hi[, 1])
    plot(
        sort(h_out$hi[, 1]),
        ylim = c(0, 1),
        pch = 19,
        xlab = "Individual (sorted by HI)",
        ylab = "Hybrid index (HI)"
    )
    # expected: columns 3 and 4 are lower/upper interval bounds
    segments(
        x0 = seq_along(ord),
        y0 = h_out$hi[ord, 3],
        x1 = seq_along(ord),
        y1 = h_out$hi[ord, 4]
    )
})

message("Fitting hierarchical genomic cline model (est_genocl) ...")
gc_out <- tryCatch(
    est_genocl(
        Gx = GlikHybrids,
        p0 = p_out$p0[, 1],
        p1 = p_out$p1[, 1],
        H = h_out$hi[, 1],
        model = "glik",
        ploidy = "diploid",
        hier = TRUE,
        n_iters = n_iters,
        p_warmup = warmup,
        n_thin = n_thin
    ),
    error = function(e) {
        message("est_genocl call failed with p_warmup/n_thin; retrying with n_iters only ...")
        est_genocl(
            Gx = GlikHybrids,
            p0 = p_out$p0[, 1],
            p1 = p_out$p1[, 1],
            H = h_out$hi[, 1],
            model = "glik",
            ploidy = "diploid",
            hier = TRUE,
            n_iters = n_iters
        )
    }
)

# save SD summaries
write_any(gc_out$SDc, file_base = paste0(prefix, "__gc_out__SDc"), outdir = text_dir)
write_any(gc_out$SDv, file_base = paste0(prefix, "__gc_out__SDv"), outdir = text_dir)

message("Applying sum-to-zero constraint (sum2zero) ...")
sz_out <- sum2zero(hmc = gc_out$gencline_hmc, transform = TRUE, ci = ci)

# ---- plots: genomic clines (raw + sum2zero) ----
gc_raw_pdf <- file.path(plots_dir, paste0(prefix, "_genclines_raw.pdf"))
save_plot_pdf(gc_raw_pdf, {
    gencline_plot(center = gc_out$center[, 1], v = gc_out$gradient[, 1], pdf = FALSE)
})

gc_sz_pdf <- file.path(plots_dir, paste0(prefix, "_genclines_sum2zero.pdf"))
save_plot_pdf(gc_sz_pdf, {
    gencline_plot(center = sz_out$center[, 1], v = sz_out$gradient[, 1], pdf = FALSE)
})

# ---- steep clines summary ----
steep_idx <- which(sz_out$gradient[, 2] > 1)
steep_n <- sum(sz_out$gradient[, 2] > 1)

write_any(steep_idx, file_base = paste0(prefix, "__sz_out__steep_gradient_indices"), outdir = text_dir)
write_any(steep_n, file_base = paste0(prefix, "__sz_out__steep_gradient_count"), outdir = text_dir)

message("Estimating interspecific ancestry (est_Q) ...")
q_out <- est_Q(
    Gx = GlikHybrids,
    p0 = p_out$p0[, 1],
    p1 = p_out$p1[, 1],
    model = "glik",
    ploidy = "diploid"
)

tri_pdf <- file.path(plots_dir, paste0(prefix, "_triangle_plot.pdf"))
save_plot_pdf(
    tri_pdf,
    {
        tri_plot(hi = q_out$hi[, 1], Q10 = q_out$Q10[, 1], pdf = FALSE, pch = 19)
    },
    width = 6.5,
    height = 6.5
)

# ---- save objects ----
sz_rdata <- file.path(rds_dir, paste0(prefix, "_sz_out.RData"))
gc_rdata <- file.path(rds_dir, paste0(prefix, "_gc_out.RData"))

save(sz_out, file = sz_rdata)
save(gc_out, file = gc_rdata)

saveRDS(
    list(p_out = p_out, h_out = h_out, gc_out = gc_out, q_out = q_out, sz_out = sz_out),
    file = file.path(rds_dir, paste0(prefix, "_all_outputs.rds"))
)

# ---- MCMC diagnostics and logs ----
message("Writing HMC draw logs and summaries ...")

write_stan_draws_wide(h_out$hi_hmc, paste0(prefix, "__hi_hmc"), text_dir)
write_stan_draws_long(h_out$hi_hmc, paste0(prefix, "__hi_hmc"), text_dir)
write_stan_summary(h_out$hi_hmc, paste0(prefix, "__hi_hmc"), text_dir)
write_stan_sampler_diagnostics(h_out$hi_hmc, paste0(prefix, "__hi_hmc"), text_dir)

write_stan_draws_wide(gc_out$gencline_hmc, paste0(prefix, "__gencline_hmc"), text_dir)
write_stan_draws_long(gc_out$gencline_hmc, paste0(prefix, "__gencline_hmc"), text_dir)
write_stan_summary(gc_out$gencline_hmc, paste0(prefix, "__gencline_hmc"), text_dir)
write_stan_sampler_diagnostics(gc_out$gencline_hmc, paste0(prefix, "__gencline_hmc"), text_dir)

# ---- write text outputs of ALL sz_out components ----
message("Writing text exports of sz_out components ...")
write_any(sz_out, file_base = paste0(prefix, "__sz_out"), outdir = text_dir)

# ---- optional: export key matrices for convenience ----
write_any(p_out$p0, file_base = paste0(prefix, "__p_out__p0"), outdir = text_dir)
write_any(p_out$p1, file_base = paste0(prefix, "__p_out__p1"), outdir = text_dir)
write_any(h_out$hi, file_base = paste0(prefix, "__h_out__hi"), outdir = text_dir)
write_any(q_out$hi, file_base = paste0(prefix, "__q_out__hi"), outdir = text_dir)
write_any(q_out$Q10, file_base = paste0(prefix, "__q_out__Q10"), outdir = text_dir)
write_any(gc_out$center, file_base = paste0(prefix, "__gc_out__center"), outdir = text_dir)
write_any(gc_out$gradient, file_base = paste0(prefix, "__gc_out__gradient"), outdir = text_dir)

message("\nDONE.")
message("Output directory: ", normalizePath(outdir, winslash = "/", mustWork = FALSE))
message("Plots: ", plots_dir)
message("Text exports: ", text_dir)
message("Saved output objects: ", gc_rdata, "  and  ", sz_rdata)
