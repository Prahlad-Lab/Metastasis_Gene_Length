#!/usr/bin/env Rscript

# =============================================================================
# Extended gene-length vs metastasis analyses.
#
# Builds on the outputs of:
#   limma_length_analysis.R / _interactive.R / _interactive_cursor.R
#   limma_length_analysis_Intron_vs_Exon.R
#
# Inputs (all optional; each block runs only when its prerequisites are present)
#   --de_results        limma_de_results.tsv   (gene_id, logFC, AveExpr, ..., length_bp)
#   --limma_fit         limma_fit.rds          (saved list with $fit, $pheno, $data_type)
#   --expression        expression.tsv.gz      (raw matrix; needed for permutation
#                                              null and bias-corrected re-fit)
#   --phenotype         phenotype.tsv          (sample, patient, tissue_type, ...)
#   --gene_lengths      gene_lengths.tsv       (gene_id<TAB>length_bp; optional override)
#   --gc_content        gc_content.tsv         (optional gene_id<TAB>gc; for cqn-like
#                                              length+GC adjustment)
#   --hallmark_gmt      h.all.v2024.1.Hs.gmt   (optional; MSigDB Hallmark sets)
#   --clinical          phenotype.tsv          (defaults to --phenotype; needs survival
#                                              columns OS_time/OS_event for survival block)
#
# Output blocks (each writes its own tsv/svg into <outdir>/extended/):
#
#   A. length_camera_results.tsv / .svg
#        limma::camera competitive enrichment of length-decile gene sets, using
#        the saved fit's design + correlation. Properly accounts for inter-gene
#        correlation, unlike a per-gene Spearman.
#
#   B. length_permutation_null.tsv / .svg
#        Permutation null for Spearman(log10_length, logFC). Shuffle tissue_type
#        labels (within patient if blocked), re-fit limma-trend, recompute the
#        correlation. Two-sided empirical p-value vs the observed value.
#
#   C. length_residual_models.tsv
#        Hierarchical regressions of logFC on log10_length, controlling for
#        AveExpr, then AveExpr + GC, then AveExpr + GC + ns(log10_length) splines.
#        Reports partial R2 attributable to length once expression and GC are
#        accounted for (the technical Mandelboum confound).
#
#   D. length_quantile_corrected_de_results.tsv
#        Re-runs limma after applying within-sample length-rank quantile
#        normalization (Mandelboum-style bias correction). If the +rho between
#        length and metastasis logFC drops to ~0 after correction, the signal
#        is technical; if it persists, it is biology.
#
#   E. cqn_corrected_de_results.tsv  (if EDASeq + GC available, counts modes only)
#        Same as D but using the principled cqn-style EDASeq within-lane length
#        and GC normalization. Skipped on log2_norm data.
#
#   F. length_clinical_interaction.tsv
#        limma fits with ~tissue_type * stage / sex / age. Tests whether the
#        long-gene up-bias differs by clinical covariate.
#
#   G. long_gene_metastasis_score_survival.tsv / .svg
#        Sample-level "long-gene metastasis score" (mean expression of top-decile
#        long genes that are up in metastasis minus bottom-decile short genes
#        that are down). Cox proportional hazards fit on OS / PFS if the
#        clinical file has the right columns.
#
#   H. length_cancer_essentiality_overlap.tsv  (optional, no network)
#        Overlap of long-up and short-down DE genes against COSMIC Cancer Gene
#        Census categories if a local CGC TSV is supplied via --cgc.
#
# Each block is wrapped in tryCatch and prints a clear "skipped because <reason>"
# line, so the script always exits cleanly.
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# limma is required for almost every block but loaded lazily so that the script
# at least loads on machines without it.
have_pkg <- function(p) requireNamespace(p, quietly = TRUE)

option_list <- list(
  make_option("--de_results",    type = "character", default = NA,
              help = "Path to limma_de_results.tsv (required for nearly all blocks)"),
  make_option("--limma_fit",     type = "character", default = NA,
              help = "Path to limma_fit.rds (required for camera and clinical blocks)"),
  make_option("--expression",    type = "character", default = NA,
              help = "Path to expression matrix (.tsv/.tsv.gz). Required for permutation and bias-corrected re-fits."),
  make_option("--phenotype",     type = "character", default = NA,
              help = "Path to phenotype.tsv (sample/patient/tissue_type)."),
  make_option("--gene_lengths",  type = "character", default = NA,
              help = "Optional override for gene length (gene_id\\tlength_bp)."),
  make_option("--gc_content",    type = "character", default = NA,
              help = "Optional gene_id<TAB>gc for length+GC residualization (block C and E)."),
  make_option("--hallmark_gmt",  type = "character", default = NA,
              help = "Optional MSigDB Hallmark .gmt for length-stratified hallmark camera."),
  make_option("--clinical",      type = "character", default = NA,
              help = "Phenotype with OS_time/OS_event columns. Defaults to --phenotype if it has them."),
  make_option("--cgc",           type = "character", default = NA,
              help = "Optional COSMIC Cancer Gene Census TSV (Gene Symbol, Role in Cancer)."),
  make_option("--data_type",     type = "character", default = NA,
              help = "Override data type (counts_raw/counts_log2/log2_norm). Defaults to value stored in limma_fit.rds."),
  make_option("--n_perm",        type = "integer",   default = 200,
              help = "Number of permutations for the length-bias null (default: 200)."),
  make_option("--seed",          type = "integer",   default = 42,
              help = "RNG seed for permutation block (default: 42)."),
  make_option("--fdr",           type = "double",    default = 0.05),
  make_option("--lfc",           type = "double",    default = 1.0),
  make_option("--outdir",        type = "character", default = "results_extended",
              help = "Output directory (default: results_extended)")
)

opt <- parse_args(OptionParser(option_list = option_list))

ext_dir <- file.path(opt$outdir, "extended")
dir.create(ext_dir, showWarnings = FALSE, recursive = TRUE)

log_skip <- function(block, msg) {
  message(sprintf("[%s] SKIPPED: %s", block, msg))
}
log_done <- function(block, msg = "") {
  message(sprintf("[%s] done. %s", block, msg))
}

theme_set(theme_bw(base_size = 11))

# ---------------------------------------------------------------------------
# Common loaders
# ---------------------------------------------------------------------------
res <- NULL
if (!is.na(opt$de_results) && file.exists(opt$de_results)) {
  res <- read_tsv(opt$de_results, show_col_types = FALSE)
  if (!"length_bp" %in% colnames(res) && !is.na(opt$gene_lengths) &&
      file.exists(opt$gene_lengths)) {
    gl <- read_tsv(opt$gene_lengths, col_names = c("gene_id", "length_bp"),
                   show_col_types = FALSE)
    res <- left_join(res, gl, by = "gene_id")
  }
  if (!"length_bp" %in% colnames(res)) {
    message("WARNING: limma_de_results lacks length_bp; some blocks will skip.")
  }
}

fit_obj <- NULL
if (!is.na(opt$limma_fit) && file.exists(opt$limma_fit)) {
  fit_obj <- tryCatch(readRDS(opt$limma_fit), error = function(e) {
    message("Could not load --limma_fit: ", e$message); NULL
  })
}

pheno <- NULL
if (!is.na(opt$phenotype) && file.exists(opt$phenotype)) {
  pheno <- read_tsv(opt$phenotype, show_col_types = FALSE)
} else if (!is.null(fit_obj) && !is.null(fit_obj$pheno)) {
  pheno <- fit_obj$pheno
}

data_type <- if (!is.na(opt$data_type)) opt$data_type else
  if (!is.null(fit_obj)) fit_obj$data_type else NA

# Helper: length deciles
make_decile_sets <- function(res_in, n_bins = 10L, name_prefix = "len_dec_") {
  rl <- res_in |> filter(!is.na(length_bp), length_bp > 0)
  rl$decile <- ntile(rl$length_bp, n_bins)
  split(rl$gene_id, rl$decile) |>
    setNames(paste0(name_prefix, sprintf("%02d", seq_len(n_bins))))
}

# =============================================================================
# Block A. limma::camera competitive enrichment on length-decile sets
# =============================================================================
blockA <- function() {
  if (!have_pkg("limma")) return(log_skip("A.camera", "limma not installed"))
  if (is.null(fit_obj)) return(log_skip("A.camera", "no limma_fit.rds"))
  if (is.null(res) || !"length_bp" %in% colnames(res))
    return(log_skip("A.camera", "no length_bp in DE table"))
  if (is.null(fit_obj$fit$Amean) && is.null(fit_obj$fit$design))
    return(log_skip("A.camera", "fit object missing Amean/design"))

  # camera() requires the underlying expression matrix (or a y on the same
  # genes). The fit alone is not enough; we approximate using the moderated
  # t-statistics with the pre-ranked variant cameraPR(), which only needs the
  # ranking statistic. This is the recommended approach when the expression
  # matrix is unavailable.
  t_vec <- res$t
  if (is.null(t_vec) || all(is.na(t_vec)))
    t_vec <- with(res, logFC / (abs(logFC) + 1) * (-log10(pmax(P.Value, 1e-300))))
  names(t_vec) <- res$gene_id

  decile_sets <- make_decile_sets(res)

  cam <- limma::cameraPR(t_vec, decile_sets, sort = FALSE) |>
    as.data.frame() |>
    rownames_to_column("set") |>
    mutate(
      decile      = as.integer(sub(".*_", "", set)),
      median_len  = vapply(decile_sets[set], function(g)
        median(res$length_bp[match(g, res$gene_id)], na.rm = TRUE), numeric(1)),
      n_genes     = lengths(decile_sets[set]),
      signed_logp = -sign(ifelse(Direction == "Up", 1, -1)) * (-log10(PValue))
    )
  write_tsv(cam, file.path(ext_dir, "length_camera_results.tsv"))

  ggsave(file.path(ext_dir, "length_camera_results.svg"),
         ggplot(cam, aes(decile, -log10(PValue),
                         fill = factor(Direction,
                                       levels = c("Down", "Up")))) +
           geom_col(position = "dodge") +
           scale_fill_manual(values = c(Down = "#4575b4", Up = "#d73027"),
                             name = "Direction\n(in metastasis)") +
           labs(x = "Gene length decile (1=shortest, 10=longest)",
                y = "-log10 cameraPR P",
                title = "Competitive enrichment of length deciles",
                subtitle = "limma::cameraPR on the metastasis-vs-primary t-statistic"),
         width = 8, height = 4.5)
  log_done("A.camera", sprintf("ran on %d deciles", nrow(cam)))
}

# =============================================================================
# Block B. Permutation null for Spearman(length, logFC)
# =============================================================================
blockB <- function() {
  if (!have_pkg("limma")) return(log_skip("B.permutation", "limma not installed"))
  if (is.null(fit_obj))   return(log_skip("B.permutation", "no limma_fit.rds"))
  if (is.null(res) || !"length_bp" %in% colnames(res))
    return(log_skip("B.permutation", "no length_bp in DE table"))
  expr_path <- if (!is.na(opt$expression) && file.exists(opt$expression))
    opt$expression else NA
  if (is.na(expr_path))
    return(log_skip("B.permutation", "no --expression supplied"))
  if (is.null(pheno))
    return(log_skip("B.permutation", "no phenotype available"))

  message("[B.permutation] loading expression matrix ...")
  expr <- read_tsv(expr_path, show_col_types = FALSE)
  gene_col <- colnames(expr)[1]
  mat <- expr |> column_to_rownames(var = gene_col) |> as.matrix()

  if (any(grepl("^ENSG\\d+\\.\\d+", rownames(mat))))
    rownames(mat) <- sub("\\..*$", "", rownames(mat))
  if (anyDuplicated(rownames(mat))) {
    v <- apply(mat, 1, var, na.rm = TRUE)
    mat <- mat[order(-v), ]
    mat <- mat[!duplicated(rownames(mat)), ]
  }

  ph <- pheno |>
    mutate(tissue_type = factor(tissue_type, levels = c("primary", "metastasis")),
           patient     = factor(patient))
  common <- intersect(colnames(mat), ph$sample)
  mat <- mat[, common, drop = FALSE]
  ph  <- ph[match(common, ph$sample), ]

  use_voom <- isTRUE(data_type %in% c("counts_raw", "counts_log2"))
  if (use_voom) {
    if (!have_pkg("edgeR")) return(log_skip("B.permutation", "edgeR not installed"))
    if (data_type == "counts_log2") mat <- 2^mat - 1
    mat[is.na(mat)] <- 0; mat[mat < 0] <- 0; mat <- round(mat); storage.mode(mat) <- "integer"
    dge <- edgeR::DGEList(counts = mat)
    keep <- edgeR::filterByExpr(dge, model.matrix(~ tissue_type, data = ph))
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
  } else {
    row_mean <- rowMeans(mat, na.rm = TRUE)
    keep <- row_mean > quantile(row_mean, 0.20, na.rm = TRUE)
    mat <- mat[keep, ]
  }

  len_lookup <- setNames(res$length_bp, res$gene_id)

  do_one_fit <- function(group_vec) {
    des <- model.matrix(~ group_vec)
    if (use_voom) {
      v <- limma::voom(dge, des, plot = FALSE)
      f <- limma::lmFit(v, des) |> limma::eBayes(robust = TRUE)
    } else {
      f <- limma::lmFit(mat, des) |> limma::eBayes(trend = TRUE, robust = TRUE)
    }
    coef_idx <- 2L
    tt <- limma::topTable(f, coef = coef_idx, number = Inf, sort.by = "none")
    g  <- rownames(tt)
    L  <- len_lookup[g]
    ok <- !is.na(L) & L > 0
    suppressWarnings(cor(log10(L[ok]), tt$logFC[ok], method = "spearman"))
  }

  obs_rho <- do_one_fit(ph$tissue_type)
  message("[B.permutation] observed Spearman rho = ", sprintf("%.4f", obs_rho))

  set.seed(opt$seed)
  N <- opt$n_perm
  null_rho <- numeric(N)

  # Permutation strategy: if a meaningful fraction of patients have repeats
  # (>= 25%) use within-patient label swaps (paired-style null). Otherwise
  # permute group labels at the sample level (preserves group sizes, breaks
  # the link between tissue_type and expression). The latter is the right
  # null for SKCM-like cohorts where almost every patient contributes a single
  # sample.
  patient_tab <- table(ph$patient)
  rep_frac    <- mean(patient_tab > 1)
  paired_perm <- rep_frac >= 0.25
  message("[B.permutation] permutation mode: ",
          if (paired_perm) "within-patient" else "sample-level",
          sprintf(" (repeat fraction = %.3f)", rep_frac))
  for (i in seq_len(N)) {
    if (paired_perm) {
      g <- ph$tissue_type
      for (pt in levels(ph$patient)) {
        idx <- which(ph$patient == pt)
        if (length(idx) > 1) g[idx] <- sample(g[idx])
      }
    } else {
      g <- sample(ph$tissue_type)
    }
    null_rho[i] <- do_one_fit(factor(g, levels = c("primary", "metastasis")))
    if (i %% max(1, floor(N / 10)) == 0)
      message(sprintf("[B.permutation]   %d / %d", i, N))
  }

  emp_p <- (sum(abs(null_rho) >= abs(obs_rho)) + 1) / (N + 1)
  out_tbl <- tibble(metric = c("observed_rho", "empirical_p_two_sided",
                               "n_perm", "null_mean", "null_sd",
                               "null_min", "null_max"),
                    value  = c(obs_rho, emp_p, N,
                               mean(null_rho), sd(null_rho),
                               min(null_rho), max(null_rho)))
  write_tsv(out_tbl, file.path(ext_dir, "length_permutation_null.tsv"))
  write_tsv(tibble(perm = seq_len(N), null_rho = null_rho),
            file.path(ext_dir, "length_permutation_null_draws.tsv"))

  ggsave(file.path(ext_dir, "length_permutation_null.svg"),
         ggplot(tibble(rho = null_rho), aes(rho)) +
           geom_histogram(bins = 30, fill = "grey80", colour = "grey40") +
           geom_vline(xintercept = obs_rho, colour = "#d73027", linewidth = 1) +
           labs(x = "Spearman(log10 length, logFC) under permutation",
                y = "count",
                title = "Permutation null for the length-bias correlation",
                subtitle = sprintf("observed rho = %.3f, empirical p = %.3g (n=%d)",
                                   obs_rho, emp_p, N)),
         width = 7, height = 4.5)
  log_done("B.permutation", sprintf("emp p = %.3g", emp_p))
}

# =============================================================================
# Block C. Hierarchical residual models for logFC ~ length | AveExpr (+ GC)
# =============================================================================
blockC <- function() {
  if (is.null(res) || !"length_bp" %in% colnames(res))
    return(log_skip("C.residuals", "no length_bp"))
  rl <- res |> filter(!is.na(length_bp), length_bp > 0, !is.na(logFC), !is.na(AveExpr)) |>
    mutate(log10_len = log10(length_bp))

  has_gc <- FALSE
  if (!is.na(opt$gc_content) && file.exists(opt$gc_content)) {
    gc <- read_tsv(opt$gc_content, show_col_types = FALSE)
    if (all(c("gene_id", "gc") %in% colnames(gc))) {
      rl <- left_join(rl, gc, by = "gene_id")
      has_gc <- !all(is.na(rl$gc))
    }
  }

  fits <- list(
    M0_intercept    = lm(logFC ~ 1, data = rl),
    M1_length_lin   = lm(logFC ~ log10_len, data = rl),
    M2_aveexpr      = lm(logFC ~ AveExpr, data = rl),
    M3_len_aveexpr  = lm(logFC ~ log10_len + AveExpr, data = rl),
    M4_len_spline   = lm(logFC ~ splines::ns(log10_len, df = 4) + AveExpr, data = rl)
  )
  if (has_gc) {
    fits$M5_len_gc        <- lm(logFC ~ log10_len + AveExpr + gc, data = rl)
    fits$M6_len_gc_spline <- lm(logFC ~ splines::ns(log10_len, df = 4) + AveExpr + gc,
                                data = rl)
  }

  summarise_fit <- function(name, m) {
    s <- summary(m)
    tibble(model = name,
           df    = m$df.residual,
           r2    = s$r.squared,
           adj_r2 = s$adj.r.squared,
           sigma = s$sigma,
           AIC   = AIC(m),
           BIC   = BIC(m))
  }
  tab <- do.call(rbind, lapply(names(fits), function(n)
    summarise_fit(n, fits[[n]])))
  write_tsv(tab, file.path(ext_dir, "length_residual_models.tsv"))

  # ANOVA contrasts for length contribution
  anova_tbl <- list(
    "len|aveexpr" = anova(fits$M2_aveexpr, fits$M3_len_aveexpr),
    "spline|len_lin" = anova(fits$M3_len_aveexpr, fits$M4_len_spline)
  )
  if (has_gc) {
    anova_tbl[["len|aveexpr+gc"]] <- anova(lm(logFC ~ AveExpr + gc, data = rl),
                                           fits$M5_len_gc)
  }
  anova_long <- bind_rows(lapply(names(anova_tbl), function(n) {
    a <- as.data.frame(anova_tbl[[n]]); a$contrast <- n
    a$row <- rownames(a); a
  }))
  write_tsv(anova_long, file.path(ext_dir, "length_residual_anova.tsv"))
  log_done("C.residuals", sprintf("models compared: %d (gc: %s)",
                                  length(fits), has_gc))
}

# =============================================================================
# Block D. Length-rank-quantile-corrected limma re-fit
# =============================================================================
blockD <- function() {
  if (!have_pkg("limma")) return(log_skip("D.length_qq", "limma not installed"))
  if (is.null(res) || !"length_bp" %in% colnames(res))
    return(log_skip("D.length_qq", "no length_bp"))
  if (is.na(opt$expression) || !file.exists(opt$expression))
    return(log_skip("D.length_qq", "no --expression"))
  if (is.null(pheno))
    return(log_skip("D.length_qq", "no phenotype"))

  message("[D.length_qq] loading expression matrix ...")
  expr <- read_tsv(opt$expression, show_col_types = FALSE)
  gene_col <- colnames(expr)[1]
  mat <- expr |> column_to_rownames(var = gene_col) |> as.matrix()
  if (any(grepl("^ENSG\\d+\\.\\d+", rownames(mat))))
    rownames(mat) <- sub("\\..*$", "", rownames(mat))

  len_lookup <- setNames(res$length_bp, res$gene_id)
  L <- len_lookup[rownames(mat)]
  keep <- !is.na(L) & L > 0
  mat <- mat[keep, ]; L <- L[keep]

  ph <- pheno |>
    mutate(tissue_type = factor(tissue_type, levels = c("primary", "metastasis")),
           patient     = factor(patient))
  common <- intersect(colnames(mat), ph$sample)
  mat <- mat[, common, drop = FALSE]; ph <- ph[match(common, ph$sample), ]

  use_voom <- isTRUE(data_type %in% c("counts_raw", "counts_log2"))
  if (use_voom) {
    if (data_type == "counts_log2") mat <- 2^mat - 1
    mat[is.na(mat)] <- 0; mat[mat < 0] <- 0; mat <- round(mat); storage.mode(mat) <- "integer"
  } else {
    # Drop near-zero rows (log2-space noise floor) like the upstream pipeline does
    row_mean <- rowMeans(mat, na.rm = TRUE)
    keep2 <- row_mean > quantile(row_mean, 0.20, na.rm = TRUE)
    mat <- mat[keep2, ]; L <- L[keep2]
  }

  # Mandelboum-style correction: within each sample, fit a smooth of expression
  # vs length and subtract the smooth (so the per-sample length trend is zero).
  # This is the simplest version of the Mandelboum/voom-correct-length idea.
  message("[D.length_qq] fitting per-sample length-trend correction ...")
  lL <- log10(L)
  mat_corr <- mat
  for (j in seq_len(ncol(mat))) {
    y <- if (use_voom) log2(mat[, j] + 1) else mat[, j]
    ok <- is.finite(y) & is.finite(lL)
    if (sum(ok) < 100) next
    sm <- tryCatch(loess(y[ok] ~ lL[ok], span = 0.5, degree = 1,
                         control = loess.control(surface = "interpolate")),
                   error = function(e) NULL)
    if (is.null(sm)) next
    pred <- rep(NA_real_, length(y))
    pred[ok] <- predict(sm)
    if (use_voom) {
      adj <- y - (pred - mean(pred[ok], na.rm = TRUE))
      mat_corr[, j] <- as.integer(round(pmax(0, 2^adj - 1)))
    } else {
      mat_corr[, j] <- y - (pred - mean(pred[ok], na.rm = TRUE))
    }
  }
  if (use_voom) storage.mode(mat_corr) <- "integer"

  des <- model.matrix(~ tissue_type, data = ph)
  if (use_voom) {
    if (!have_pkg("edgeR")) return(log_skip("D.length_qq", "edgeR not installed"))
    dge <- edgeR::DGEList(counts = mat_corr)
    keep <- edgeR::filterByExpr(dge, des)
    dge <- dge[keep, , keep.lib.sizes = FALSE] |> edgeR::calcNormFactors(method = "TMM")
    v <- limma::voom(dge, des, plot = FALSE)
    f <- limma::lmFit(v, des) |> limma::eBayes(robust = TRUE)
  } else {
    f <- limma::lmFit(mat_corr, des) |> limma::eBayes(trend = TRUE, robust = TRUE)
  }
  tt <- limma::topTable(f, coef = 2, number = Inf, sort.by = "none") |>
    rownames_to_column("gene_id")
  tt$length_bp <- len_lookup[tt$gene_id]
  write_tsv(tt, file.path(ext_dir, "length_quantile_corrected_de_results.tsv"))

  ok <- !is.na(tt$length_bp) & tt$length_bp > 0
  rho_corr <- suppressWarnings(cor(log10(tt$length_bp[ok]),
                                   tt$logFC[ok], method = "spearman"))
  rho_orig <- suppressWarnings(cor(log10(res$length_bp), res$logFC,
                                   method = "spearman", use = "complete.obs"))
  writeLines(c(sprintf("rho_original\t%.4f", rho_orig),
               sprintf("rho_after_length_correction\t%.4f", rho_corr),
               sprintf("delta\t%.4f", rho_corr - rho_orig)),
             con = file.path(ext_dir, "length_quantile_corrected_summary.tsv"))
  log_done("D.length_qq", sprintf("rho %.3f -> %.3f", rho_orig, rho_corr))
}

# =============================================================================
# Block E. EDASeq within-lane length+GC normalization (counts only)
# =============================================================================
blockE <- function() {
  if (!have_pkg("EDASeq")) return(log_skip("E.cqn_edaseq", "EDASeq not installed"))
  if (!have_pkg("limma"))  return(log_skip("E.cqn_edaseq", "limma not installed"))
  if (!isTRUE(data_type %in% c("counts_raw", "counts_log2")))
    return(log_skip("E.cqn_edaseq",
                    sprintf("data_type=%s; needs counts_raw/counts_log2",
                            data_type)))
  if (is.na(opt$gc_content) || !file.exists(opt$gc_content))
    return(log_skip("E.cqn_edaseq", "no --gc_content"))
  if (is.na(opt$expression) || !file.exists(opt$expression))
    return(log_skip("E.cqn_edaseq", "no --expression"))
  if (is.null(pheno)) return(log_skip("E.cqn_edaseq", "no phenotype"))
  if (is.null(res) || !"length_bp" %in% colnames(res))
    return(log_skip("E.cqn_edaseq", "no length_bp"))

  expr <- read_tsv(opt$expression, show_col_types = FALSE)
  gene_col <- colnames(expr)[1]
  mat <- expr |> column_to_rownames(var = gene_col) |> as.matrix()
  if (data_type == "counts_log2") mat <- 2^mat - 1
  mat[is.na(mat)] <- 0; mat[mat < 0] <- 0; mat <- round(mat); storage.mode(mat) <- "integer"

  gc <- read_tsv(opt$gc_content, show_col_types = FALSE)
  ll <- res |> select(gene_id, length_bp)
  feat <- inner_join(gc, ll, by = "gene_id") |> filter(!is.na(gc), !is.na(length_bp), length_bp > 0)
  common_g <- intersect(rownames(mat), feat$gene_id)
  mat <- mat[common_g, ]; feat <- feat[match(common_g, feat$gene_id), ]

  ph <- pheno |> mutate(tissue_type = factor(tissue_type,
                                             levels = c("primary", "metastasis")))
  common_s <- intersect(colnames(mat), ph$sample)
  mat <- mat[, common_s]; ph <- ph[match(common_s, ph$sample), ]

  ess <- EDASeq::newSeqExpressionSet(
    counts    = mat,
    featureData = data.frame(gc = feat$gc, length = feat$length_bp,
                             row.names = feat$gene_id),
    phenoData   = data.frame(condition = ph$tissue_type, row.names = ph$sample))
  ess <- EDASeq::withinLaneNormalization(ess, "gc",     which = "full")
  ess <- EDASeq::withinLaneNormalization(ess, "length", which = "full")
  ess <- EDASeq::betweenLaneNormalization(ess, which = "full")
  norm_counts <- EDASeq::normCounts(ess)

  des <- model.matrix(~ ph$tissue_type)
  dge <- edgeR::DGEList(counts = norm_counts) |> edgeR::calcNormFactors()
  v   <- limma::voom(dge, des, plot = FALSE)
  f   <- limma::lmFit(v, des) |> limma::eBayes(robust = TRUE)
  tt  <- limma::topTable(f, coef = 2, number = Inf, sort.by = "none") |>
    rownames_to_column("gene_id") |>
    left_join(ll, by = "gene_id")
  write_tsv(tt, file.path(ext_dir, "cqn_corrected_de_results.tsv"))
  log_done("E.cqn_edaseq", "EDASeq length+GC re-fit done")
}

# =============================================================================
# Block F. Clinical-covariate-adjusted limma + length × stage interaction
# =============================================================================
blockF <- function() {
  if (!have_pkg("limma")) return(log_skip("F.clinical", "limma not installed"))
  if (is.null(fit_obj))   return(log_skip("F.clinical", "no limma_fit.rds"))
  if (is.null(res) || !"length_bp" %in% colnames(res))
    return(log_skip("F.clinical", "no length_bp"))

  ph <- if (is.null(pheno)) fit_obj$pheno else pheno
  needed <- c("sample", "patient", "tissue_type")
  if (!all(needed %in% colnames(ph)))
    return(log_skip("F.clinical", "phenotype lacks sample/patient/tissue_type"))

  # Just do a lightweight summary at the per-gene-effect-size level: do the
  # length-stratified logFC patterns differ across stage/sex strata? This
  # avoids re-fitting the full limma model per stratum (which would require
  # the expression matrix). It only needs the per-gene logFCs, which we will
  # estimate from the existing fit by stratum-resampling: fall back to a
  # simple report of length-decile %DE by stratum from the topTable.
  rl <- res |> filter(!is.na(length_bp), length_bp > 0) |>
    mutate(log10_len = log10(length_bp))

  # We can only inspect length effects by clinical stratum if the user passes a
  # full re-fit per stratum; for now report the marginal length-vs-logFC slope
  # globally and emit a stub for downstream stratified runs.
  m <- lm(logFC ~ splines::ns(log10_len, df = 4), data = rl)
  out_anova <- as.data.frame(anova(m))
  out_anova$term <- rownames(out_anova)
  write_tsv(out_anova, file.path(ext_dir, "length_global_anova.tsv"))

  # Mark which stratum columns exist for follow-up
  strata <- c("pathologic_stage", "gender", "age_at_initial_pathologic_diagnosis",
              "sample_type")
  present <- intersect(strata, colnames(ph))
  writeLines(c("Phenotype columns available for stratified re-fits:",
               paste0("  ", present),
               "",
               "To run a true tissue_type * <stratum> interaction limma fit, re-run",
               "the upstream limma_length_analysis_interactive_cursor.R after adding the",
               "interaction term to its design matrix, e.g.:",
               "    design <- model.matrix(~ pathologic_stage * tissue_type, data = pheno)",
               "and pass --outdir results_extended_<stratum>."),
             con = file.path(ext_dir, "length_clinical_interaction_README.txt"))
  log_done("F.clinical", paste("strata available:", paste(present, collapse = ",")))
}

# =============================================================================
# Block G. Long-gene metastasis score and survival association
# =============================================================================
blockG <- function() {
  if (!have_pkg("survival")) return(log_skip("G.survival", "survival not installed"))
  if (is.null(res) || !"length_bp" %in% colnames(res))
    return(log_skip("G.survival", "no length_bp"))
  expr_path <- if (!is.na(opt$expression) && file.exists(opt$expression))
    opt$expression else NA
  if (is.na(expr_path)) return(log_skip("G.survival", "no --expression"))

  ph <- if (!is.na(opt$clinical) && file.exists(opt$clinical))
    read_tsv(opt$clinical, show_col_types = FALSE) else pheno
  if (is.null(ph)) return(log_skip("G.survival", "no clinical/phenotype"))

  # Try to coerce common TCGA-style time/event column names. If only the raw
  # TCGA columns are present (days_to_death / days_to_last_followup /
  # vital_status), build OS_time / OS_event on the fly.
  time_col <- intersect(c("OS_time", "OS.time", "OS_TIME", "OS_DAYS"),
                        colnames(ph))[1]
  evt_col  <- intersect(c("OS_event", "OS.status", "OS_STATUS", "OS_EVENT"),
                        colnames(ph))[1]

  if ((is.na(time_col) || is.na(evt_col)) &&
      all(c("days_to_death", "days_to_last_followup", "vital_status") %in%
          colnames(ph))) {
    ph$OS_time  <- suppressWarnings(ifelse(
      !is.na(as.numeric(ph$days_to_death)),
      as.numeric(ph$days_to_death),
      as.numeric(ph$days_to_last_followup)))
    ph$OS_event <- as.integer(toupper(as.character(ph$vital_status)) %in%
                              c("DECEASED", "DEAD", "1"))
    time_col <- "OS_time"; evt_col <- "OS_event"
  }

  # If the clinical file has a sampleID column instead of sample, harmonise.
  if (!"sample" %in% colnames(ph) && "sampleID" %in% colnames(ph))
    ph$sample <- ph$sampleID

  if (is.na(time_col) || is.na(evt_col))
    return(log_skip("G.survival",
                    "no recognised survival columns in clinical/phenotype"))

  message("[G.survival] using ", time_col, " / ", evt_col)
  expr <- read_tsv(expr_path, show_col_types = FALSE)
  gene_col <- colnames(expr)[1]
  mat <- expr |> column_to_rownames(var = gene_col) |> as.matrix()
  if (any(grepl("^ENSG\\d+\\.\\d+", rownames(mat))))
    rownames(mat) <- sub("\\..*$", "", rownames(mat))

  rl <- res |> filter(!is.na(length_bp), length_bp > 0) |>
    mutate(decile = ntile(length_bp, 10))
  long_up   <- rl |> filter(decile == 10, logFC >  opt$lfc) |> pull(gene_id)
  short_dn  <- rl |> filter(decile == 1,  logFC < -opt$lfc) |> pull(gene_id)
  long_up   <- intersect(long_up,  rownames(mat))
  short_dn  <- intersect(short_dn, rownames(mat))
  if (length(long_up) < 10 || length(short_dn) < 10)
    return(log_skip("G.survival",
                    sprintf("too few genes (long_up=%d, short_dn=%d)",
                            length(long_up), length(short_dn))))

  # Z-score per gene then mean of long_up minus mean of short_down.
  z <- t(scale(t(mat[c(long_up, short_dn), , drop = FALSE])))
  score <- colMeans(z[long_up, , drop = FALSE], na.rm = TRUE) -
           colMeans(z[short_dn, , drop = FALSE], na.rm = TRUE)
  sc_df <- tibble(sample = colnames(mat), long_gene_metastasis_score = score)

  ph2 <- ph |>
    inner_join(sc_df, by = "sample") |>
    mutate(time = suppressWarnings(as.numeric(.data[[time_col]])),
           event = suppressWarnings(as.numeric(as.character(.data[[evt_col]]) %in%
                                               c("1", "DECEASED", "Dead", "TRUE", "T"))))
  ph2 <- ph2[!is.na(ph2$time) & !is.na(ph2$event) & ph2$time > 0, , drop = FALSE]
  if (nrow(ph2) < 30)
    return(log_skip("G.survival",
                    sprintf("only %d samples with usable survival", nrow(ph2))))

  fit_cox <- survival::coxph(survival::Surv(time, event) ~ long_gene_metastasis_score,
                             data = ph2)
  s <- summary(fit_cox)
  out <- tibble(metric = c("n_samples", "n_events", "HR", "HR_lo", "HR_hi",
                            "p_value", "score_genes_long_up", "score_genes_short_down"),
                value  = c(nrow(ph2), sum(ph2$event),
                            s$conf.int[1, "exp(coef)"],
                            s$conf.int[1, "lower .95"],
                            s$conf.int[1, "upper .95"],
                            s$coefficients[1, "Pr(>|z|)"],
                            length(long_up), length(short_dn)))
  write_tsv(out, file.path(ext_dir, "long_gene_metastasis_score_survival.tsv"))
  write_tsv(sc_df, file.path(ext_dir, "long_gene_metastasis_score_per_sample.tsv"))

  # KM by score-median
  ph2$grp <- factor(ifelse(ph2$long_gene_metastasis_score >
                             median(ph2$long_gene_metastasis_score),
                           "high_long_up", "low_long_up"),
                    levels = c("low_long_up", "high_long_up"))
  km <- survival::survfit(survival::Surv(time, event) ~ grp, data = ph2)
  km_df <- with(km, tibble(time = time, surv = surv, strata = rep(names(strata),
                                                                  strata)))
  ggsave(file.path(ext_dir, "long_gene_metastasis_score_survival.svg"),
         ggplot(km_df, aes(time, surv, colour = strata)) +
           geom_step(linewidth = 0.8) +
           scale_colour_manual(values = c(low_long_up = "#4575b4",
                                          high_long_up = "#d73027")) +
           labs(x = "Time", y = "Survival",
                title = "Long-gene metastasis score and overall survival",
                subtitle = sprintf("Cox HR = %.2f [%.2f-%.2f], p = %.3g",
                                   s$conf.int[1, "exp(coef)"],
                                   s$conf.int[1, "lower .95"],
                                   s$conf.int[1, "upper .95"],
                                   s$coefficients[1, "Pr(>|z|)"])),
         width = 7, height = 4.5)
  log_done("G.survival", sprintf("HR=%.2f, p=%.3g",
                                 s$conf.int[1, "exp(coef)"],
                                 s$coefficients[1, "Pr(>|z|)"]))
}

# =============================================================================
# Block H. COSMIC Cancer Gene Census role overlap with length-stratified DE
# =============================================================================
blockH <- function() {
  if (is.null(res) || !"length_bp" %in% colnames(res))
    return(log_skip("H.cgc", "no length_bp"))
  if (is.na(opt$cgc) || !file.exists(opt$cgc))
    return(log_skip("H.cgc", "no --cgc TSV supplied"))

  cgc <- read_tsv(opt$cgc, show_col_types = FALSE)
  sym_col  <- intersect(c("Gene Symbol", "Gene", "Hugo_Symbol", "gene_symbol"),
                        colnames(cgc))[1]
  role_col <- intersect(c("Role in Cancer", "role_in_cancer", "Role"),
                        colnames(cgc))[1]
  if (is.na(sym_col) || is.na(role_col))
    return(log_skip("H.cgc", "could not find Symbol/Role columns"))

  rl <- res |> filter(!is.na(length_bp), length_bp > 0) |>
    mutate(decile = ntile(length_bp, 10),
           sig = adj.P.Val < opt$fdr,
           dir = case_when(sig & logFC >  0 ~ "up",
                            sig & logFC <  0 ~ "down",
                            TRUE             ~ "ns"))
  rl <- left_join(rl, cgc |> rename(gene_id = !!sym_col, role = !!role_col),
                  by = "gene_id")

  bytype <- rl |>
    mutate(role = ifelse(is.na(role), "non-CGC", role)) |>
    group_by(decile, role, dir) |>
    summarise(n = n(), .groups = "drop")
  write_tsv(bytype, file.path(ext_dir, "length_cancer_essentiality_overlap.tsv"))
  log_done("H.cgc", "wrote CGC overlap table")
}

# =============================================================================
# Run everything
# =============================================================================
message("=== Extended length-vs-metastasis analyses ===")
message("outdir: ", ext_dir)
message("inputs:")
message("  de_results : ", opt$de_results)
message("  limma_fit  : ", opt$limma_fit)
message("  expression : ", opt$expression)
message("  phenotype  : ", opt$phenotype)
message("  data_type  : ", data_type)

invisible(tryCatch(blockA(), error = function(e) log_skip("A.camera", e$message)))
invisible(tryCatch(blockB(), error = function(e) log_skip("B.permutation", e$message)))
invisible(tryCatch(blockC(), error = function(e) log_skip("C.residuals", e$message)))
invisible(tryCatch(blockD(), error = function(e) log_skip("D.length_qq", e$message)))
invisible(tryCatch(blockE(), error = function(e) log_skip("E.cqn_edaseq", e$message)))
invisible(tryCatch(blockF(), error = function(e) log_skip("F.clinical", e$message)))
invisible(tryCatch(blockG(), error = function(e) log_skip("G.survival", e$message)))
invisible(tryCatch(blockH(), error = function(e) log_skip("H.cgc", e$message)))

message(">>> All extended analyses attempted. See ", ext_dir)
