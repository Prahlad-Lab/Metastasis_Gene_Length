#!/usr/bin/env Rscript

# =============================================================================
# limma DE + gene-length bias analysis, matrix-input (Xena) version.
#
# Data types supported:
#   counts_raw  -> round and run voomWithQualityWeights
#   counts_log2 -> reverse transform (2^x - 1), round, then voom
#   log2_norm   -> use limma-trend directly (voom not appropriate)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(limma)
  library(edgeR)
  library(goseq)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
})

# ---------------------------------------------------------------------------
# Parameter Setup (Handles both RStudio Interactive and CLI Execution)
# ---------------------------------------------------------------------------
if (interactive()) {
  # --- RSTUDIO INTERACTIVE MODE ---
  # Manually set your parameters here when running line-by-line in RStudio
  opt <- list(
    expression   = "/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/results/xena_download/expression.tsv.gz", # UPDATE THIS PATH
    phenotype    = "/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/results/xena_download/phenotype.tsv",    # UPDATE THIS PATH
    data_type    = "log2_norm",                         # "counts_raw" | "counts_log2" | "log2_norm"
    gene_lengths = "/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/gene_lengths_interactive_test.tsv",                                   # Set to path/to/gene_lengths.tsv if you have one, or keep NA
    genome       = "hg38",
    gene_id_type = "geneSymbol",
    fdr          = 0.05,
    lfc          = 1.0,
    outdir       = "/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/results_interactive_cursor"
  )
  message(">>> Running interactively in RStudio. Using hardcoded parameters.")
} else {
  # --- COMMAND LINE MODE ---
  option_list <- list(
    make_option("--expression",   type = "character"),
    make_option("--phenotype",    type = "character"),
    make_option("--data_type",    type = "character",
                help = "counts_raw | counts_log2 | log2_norm"),
    make_option("--gene_lengths", type = "character", default = NA),
    make_option("--genome",       type = "character", default = "hg38"),
    make_option("--gene_id_type", type = "character", default = "geneSymbol"),
    make_option("--fdr",          type = "double",    default = 0.05),
    make_option("--lfc",          type = "double",    default = 1.0),
    make_option("--outdir",       type = "character", default = "results")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load expression + phenotype
# =============================================================================
message(">>> Loading expression: ", opt$expression)
expr <- read_tsv(opt$expression, show_col_types = FALSE)
gene_col <- colnames(expr)[1]
mat <- expr |> column_to_rownames(var = gene_col) |> as.matrix()
n_genes_input <- nrow(mat)

# Strip Ensembl version suffixes if present (e.g. ENSG00000123456.7 -> ENSG00000123456)
if (any(grepl("^ENSG\\d+\\.\\d+", rownames(mat)))) {
  rownames(mat) <- sub("\\..*$", "", rownames(mat))
  # Collapse duplicates by taking the max-variance row
  if (anyDuplicated(rownames(mat))) {
    v <- apply(mat, 1, var, na.rm = TRUE)
    mat <- mat[order(-v), ]
    mat <- mat[!duplicated(rownames(mat)), ]
  }
}

pheno <- read_tsv(opt$phenotype, show_col_types = FALSE) |>
  mutate(
    tissue_type = factor(tissue_type, levels = c("primary", "metastasis")),
    patient     = factor(patient)
  )

# Align samples
common <- intersect(colnames(mat), pheno$sample)
mat   <- mat[, common, drop = FALSE]
pheno <- pheno[match(common, pheno$sample), ]
stopifnot(all(colnames(mat) == pheno$sample))
message("    ", nrow(mat), " genes x ", ncol(mat), " samples")
message("    primary: ", sum(pheno$tissue_type == "primary"),
        ",  metastasis: ", sum(pheno$tissue_type == "metastasis"))

# Repeated measures support: use duplicateCorrelation whenever patient IDs repeat.
patient_tab <- table(pheno$patient)
has_repeats <- any(patient_tab > 1)
has_both_groups <- all(c("primary", "metastasis") %in% as.character(unique(pheno$tissue_type)))
use_block <- has_repeats && has_both_groups
message(">>> Design: ",
        if (use_block) {
          sprintf("repeated-measures (~tissue_type; block=patient, %d repeated patients)",
                  sum(patient_tab > 1))
        } else {
          "unpaired (~tissue_type)"
        })

# =============================================================================
# 2. Branch on data type
# =============================================================================
use_voom <- opt$data_type %in% c("counts_raw", "counts_log2")
n_genes_after_filter <- NA_integer_

if (use_voom) {
  message(">>> Preparing counts for voom")
  if (opt$data_type == "counts_log2") {
    mat <- 2 ^ mat - 1
  }
  # NA handling + clip to non-negative integers
  mat[is.na(mat)] <- 0
  mat[mat < 0]    <- 0
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  
  dge <- DGEList(counts = mat)
  design_filter <- model.matrix(~ tissue_type, data = pheno)
  keep <- filterByExpr(dge, design_filter)
  message("    filterByExpr kept ", sum(keep), " / ", length(keep), " genes")
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  n_genes_after_filter <- nrow(dge)
  dge  <- calcNormFactors(dge, method = "TMM")
}

# =============================================================================
# 3. Design and fit
# =============================================================================
design <- model.matrix(~ tissue_type, data = pheno)
block <- if (use_block) pheno$patient else NULL
dupcor <- NULL

if (use_voom) {
  svg(file.path(opt$outdir, "voom_mean_variance.svg"),
      width = 800 / 96, height = 600 / 96)
  if (use_block) {
    v0 <- voom(dge, design, plot = FALSE)
    dupcor <- duplicateCorrelation(v0, design, block = block)
    message("    duplicateCorrelation consensus = ",
            sprintf("%.4f", dupcor$consensus))
    v <- voomWithQualityWeights(dge, design, block = block,
                                correlation = dupcor$consensus, plot = TRUE)
  } else {
    v <- voomWithQualityWeights(dge, design, plot = TRUE)
  }
  dev.off()
  if (use_block) {
    dupcor <- duplicateCorrelation(v, design, block = block)
    fit <- lmFit(v, design, block = block, correlation = dupcor$consensus)
  } else {
    fit <- lmFit(v, design)
  }
} else {
  message(">>> Using limma-trend on log2-normalised expression")
  # Drop near-zero rows (log2-space noise floor)
  row_mean <- rowMeans(mat, na.rm = TRUE)
  keep <- row_mean > quantile(row_mean, 0.20, na.rm = TRUE)
  mat  <- mat[keep, ]
  n_genes_after_filter <- nrow(mat)
  message("    kept ", nrow(mat), " genes after noise-floor filter")
  v   <- mat   # for downstream uniformity
  if (use_block) {
    dupcor <- duplicateCorrelation(mat, design, block = block)
    message("    duplicateCorrelation consensus = ",
            sprintf("%.4f", dupcor$consensus))
    fit <- lmFit(mat, design, block = block, correlation = dupcor$consensus)
  } else {
    fit <- lmFit(mat, design)
  }
}

fit <- eBayes(fit, trend = !use_voom, robust = TRUE)
fit_treat <- treat(fit, lfc = opt$lfc)

coef_name <- "tissue_typemetastasis"
stopifnot(coef_name %in% colnames(fit$coefficients))

res <- topTable(fit_treat, coef = coef_name, number = Inf, sort.by = "none") |>
  rownames_to_column("gene_id")

# =============================================================================
# 4. Gene lengths
# =============================================================================
if (!is.na(opt$gene_lengths)) {
  gene_len <- read_tsv(opt$gene_lengths,
                       col_names = c("gene_id", "length_bp"),
                       show_col_types = FALSE)
  message(">>> Using user-supplied lengths (n=", nrow(gene_len), ")")
} else {
  message(">>> Using goseq::getlength(", opt$genome, ", ", opt$gene_id_type, ")")
  len_vec <- tryCatch(
    suppressWarnings(getlength(res$gene_id, opt$genome, opt$gene_id_type)),
    error = function(e) {
      message("    goseq::getlength failed: ", e$message)
      NULL
    }
  )
  
  if (is.null(len_vec) || all(is.na(len_vec))) {
    message(">>> Falling back to biomaRt to retrieve gene lengths")
    len_vec <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("biomaRt not installed")
      mart <- biomaRt::useMart("ensembl",
                               dataset = "hsapiens_gene_ensembl",
                               host    = "https://ensembl.org")
      bm <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "transcript_length"),
        filters    = "hgnc_symbol",
        values     = res$gene_id,
        mart       = mart
      )
      # Use the longest transcript length per gene
      bm_agg  <- bm |>
        group_by(hgnc_symbol) |>
        summarise(length_bp = max(transcript_length, na.rm = TRUE),
                  .groups = "drop")
      len_map <- setNames(bm_agg$length_bp, bm_agg$hgnc_symbol)
      len_map[res$gene_id]
    }, error = function(e) {
      message("    biomaRt fallback failed: ", e$message)
      message("    Gene lengths unavailable; length-stratified analyses will be skipped.")
      rep(NA_real_, nrow(res))
    })
  }
  gene_len <- tibble(gene_id = res$gene_id, length_bp = len_vec)
}
res <- res |> left_join(gene_len, by = "gene_id")

write_tsv(res, file.path(opt$outdir, "limma_de_results.tsv"))
saveRDS(list(fit = fit, pheno = pheno, res = res, data_type = opt$data_type),
        file.path(opt$outdir, "limma_fit.rds"))

# =============================================================================
# 5a. Basic DE plots (always generated, no gene lengths required)
# =============================================================================
message(">>> Generating basic DE plots")
theme_set(theme_bw(base_size = 11))

res_de <- res |>
  mutate(
    sig       = adj.P.Val < opt$fdr,
    direction = case_when(sig & logFC >  0 ~ "up",
                          sig & logFC <  0 ~ "down",
                          TRUE                   ~ "ns")
  )

ggsave(file.path(opt$outdir, "volcano.svg"),
       ggplot(res_de, aes(logFC, -log10(adj.P.Val), colour = direction)) +
         geom_point(size = 0.6, alpha = 0.5) +
         geom_vline(xintercept = c(-opt$lfc, opt$lfc), linetype = 2) +
         geom_hline(yintercept = -log10(opt$fdr),        linetype = 2) +
         scale_colour_manual(values = c(up = "#d73027", down = "#4575b4", ns = "grey70")) +
         labs(title = "Volcano plot", x = "log2 FC (metastasis vs primary)",
              y = "-log10 adj.P",
              subtitle = sprintf("Up: %d  Down: %d  (treat test; FDR<%.2f, lfc=%.1f)",
                                 sum(res_de$direction == "up"),
                                 sum(res_de$direction == "down"),
                                 opt$fdr, opt$lfc)),
       width = 7, height = 5)

ggsave(file.path(opt$outdir, "ma_plot.svg"),
       ggplot(res_de, aes(AveExpr, logFC, colour = direction)) +
         geom_point(size = 0.5, alpha = 0.4) +
         geom_hline(yintercept = c(-opt$lfc, 0, opt$lfc),
                    linetype = c(2, 1, 2), colour = c("grey40", "black", "grey40")) +
         scale_colour_manual(values = c(up = "#d73027", down = "#4575b4", ns = "grey70")) +
         labs(title = "MA plot", x = "Average expression (log2)",
              y = "log2 FC (metastasis vs primary)"),
       width = 7, height = 5)

# =============================================================================
# 5b. Length-stratified summaries + plots
# =============================================================================
lengths_available <- sum(!is.na(res$length_bp)) >= 10L
n_lengths_nonmissing <- sum(!is.na(res$length_bp))
n_complete_pairs <- sum(tapply(pheno$tissue_type, pheno$patient,
                               function(x) length(unique(x)) == 2))

diag_tbl <- tibble(
  metric = c(
    "data_type", "use_voom", "design_mode", "n_samples", "n_primary", "n_metastasis",
    "n_patients", "n_repeated_patients", "n_complete_pairs", "use_block",
    "duplicate_correlation_consensus", "n_genes_input", "n_genes_after_filter",
    "n_genes_tested", "fdr_threshold", "lfc_threshold_treat",
    "n_lengths_nonmissing", "pct_lengths_nonmissing", "lengths_available"
  ),
  value = c(
    as.character(opt$data_type),
    as.character(use_voom),
    ifelse(use_block, "repeated_measures_blocked", "unpaired"),
    as.character(ncol(mat)),
    as.character(sum(pheno$tissue_type == "primary")),
    as.character(sum(pheno$tissue_type == "metastasis")),
    as.character(nlevels(pheno$patient)),
    as.character(sum(patient_tab > 1)),
    as.character(n_complete_pairs),
    as.character(use_block),
    ifelse(is.null(dupcor), NA_character_, sprintf("%.6f", dupcor$consensus)),
    as.character(n_genes_input),
    as.character(n_genes_after_filter),
    as.character(nrow(res)),
    as.character(opt$fdr),
    as.character(opt$lfc),
    as.character(n_lengths_nonmissing),
    sprintf("%.2f", 100 * n_lengths_nonmissing / max(1, nrow(res))),
    as.character(lengths_available)
  )
)
write_tsv(diag_tbl, file.path(opt$outdir, "analysis_diagnostics.tsv"))
writeLines(
  c(
    "Run diagnostics",
    paste0("design_mode: ", ifelse(use_block, "repeated_measures_blocked", "unpaired")),
    paste0("samples: ", ncol(mat), " (primary=", sum(pheno$tissue_type == "primary"),
           ", metastasis=", sum(pheno$tissue_type == "metastasis"), ")"),
    paste0("patients: ", nlevels(pheno$patient), " (repeated=", sum(patient_tab > 1),
           ", complete_pairs=", n_complete_pairs, ")"),
    paste0("duplicate_correlation_consensus: ",
           ifelse(is.null(dupcor), "NA", sprintf("%.6f", dupcor$consensus))),
    paste0("genes: input=", n_genes_input, ", after_filter=", n_genes_after_filter,
           ", tested=", nrow(res)),
    paste0("lengths: nonmissing=", n_lengths_nonmissing, " / ", nrow(res),
           " (", sprintf("%.2f", 100 * n_lengths_nonmissing / max(1, nrow(res))), "%)"),
    paste0("thresholds: FDR<", opt$fdr, ", treat_lfc=", opt$lfc)
  ),
  con = file.path(opt$outdir, "analysis_diagnostics.txt")
)

if (!lengths_available) {
  message(">>> Skipping length-stratified analysis: insufficient gene length data")
} else {
  
  res_len <- res |>
    filter(!is.na(length_bp), length_bp > 0) |>
    mutate(
      log10_len = log10(length_bp),
      decile    = ntile(length_bp, 10),
      sig       = adj.P.Val < opt$fdr,
      direction = case_when(sig & logFC >  0 ~ "up",
                            sig & logFC <  0 ~ "down",
                            TRUE                   ~ "ns")
    )
  
  sp <- suppressWarnings(cor.test(res_len$log10_len, res_len$logFC,
                                  method = "spearman"))
  cat(sprintf(">>> Spearman cor(log10 length, logFC) = %.3f, p = %.2g\n",
              sp$estimate, sp$p.value))
  writeLines(c(sprintf("spearman_rho\t%.4f", sp$estimate),
               sprintf("p_value\t%.3e",      sp$p.value),
               sprintf("n\t%d",               nrow(res_len))),
             con = file.path(opt$outdir, "length_logfc_correlation.tsv"))
  
  decile_summary <- res_len |>
    group_by(decile) |>
    summarise(n             = n(),
              median_length = median(length_bp),
              median_logFC  = median(logFC),
              pct_sig       = 100 * mean(sig),
              pct_up        = 100 * mean(direction == "up"),
              pct_down      = 100 * mean(direction == "down"),
              .groups       = "drop")
  write_tsv(decile_summary, file.path(opt$outdir, "length_decile_summary.tsv"))

  # Adjusted model: non-linear length effect while controlling for mean expression.
  len_model <- lm(logFC ~ splines::ns(log10_len, df = 3) + AveExpr, data = res_len)
  len_terms <- rownames(coef(summary(len_model)))
  len_term_idx <- grepl("ns\\(log10_len", len_terms)
  len_coef <- coef(summary(len_model))[len_term_idx, , drop = FALSE] |>
    as.data.frame() |>
    rownames_to_column("term")
  len_drop <- drop1(len_model, test = "F") |>
    as.data.frame() |>
    rownames_to_column("term")
  write_tsv(len_coef, file.path(opt$outdir, "length_adjusted_model_coefficients.tsv"))
  write_tsv(len_drop, file.path(opt$outdir, "length_adjusted_model_drop1.tsv"))

  pred_df <- tibble(
    log10_len = seq(min(res_len$log10_len), max(res_len$log10_len), length.out = 200),
    AveExpr = median(res_len$AveExpr, na.rm = TRUE)
  )
  pred_df$pred_logFC <- predict(len_model, newdata = pred_df)
  
  ggsave(file.path(opt$outdir, "length_vs_logfc.svg"),
         ggplot(res_len, aes(log10_len, logFC)) +
           geom_hex(bins = 60) +
           geom_hline(yintercept = 0, linetype = 2) +
          geom_smooth(method = "loess", se = FALSE, colour = "red") +
          geom_line(data = pred_df, aes(x = log10_len, y = pred_logFC),
                    inherit.aes = FALSE, colour = "black", linewidth = 1.0) +
           scale_fill_viridis_c(trans = "log10") +
           labs(x = "log10 gene length (bp)",
                y = "log2 FC (metastasis vs primary)",
                title = "Length vs fold change",
                subtitle = sprintf("Spearman rho = %.3f, p = %.2g",
                                   sp$estimate, sp$p.value)),
         width = 7, height = 5)
  
  ggsave(file.path(opt$outdir, "length_decile_de_pct.svg"),
         decile_summary |>
           select(decile, pct_up, pct_down) |>
           pivot_longer(c(pct_up, pct_down), names_to = "dir", values_to = "pct") |>
           mutate(dir = recode(dir, pct_up = "Up in met", pct_down = "Down in met")) |>
           ggplot(aes(factor(decile), pct, fill = dir)) +
           geom_col(position = "dodge") +
           scale_fill_manual(values = c("Up in met"   = "#d73027",
                                        "Down in met" = "#4575b4")) +
           labs(x = "Gene length decile (1=shortest, 10=longest)",
                y = sprintf("%% DE (FDR<%.2f, |logFC|>%.1f)", opt$fdr, opt$lfc),
                fill = NULL, title = "DE enrichment by length decile"),
         width = 8, height = 4.5)
  
  ggsave(file.path(opt$outdir, "length_density_by_direction.svg"),
         ggplot(res_len, aes(log10_len, colour = direction, fill = direction)) +
           geom_density(alpha = 0.15) +
           scale_colour_manual(values = c(up = "#d73027", down = "#4575b4", ns = "grey50")) +
           scale_fill_manual  (values = c(up = "#d73027", down = "#4575b4", ns = "grey50")) +
           labs(x = "log10 gene length (bp)", y = "density",
                title = "Gene length distribution by DE direction"),
         width = 7, height = 4)
  
} # end lengths_available

# =============================================================================
# 6. goseq: length-bias-corrected GO/KEGG
# =============================================================================
if (!lengths_available) {
  message(">>> Skipping goseq enrichment: insufficient gene length data")
} else {
  
  message(">>> goseq enrichment")
  gene_vec <- as.integer(res$adj.P.Val < opt$fdr)
  names(gene_vec) <- res$gene_id
  gene_vec[is.na(gene_vec)] <- 0
  bias <- setNames(res$length_bp, res$gene_id)
  bias[is.na(bias)] <- median(bias, na.rm = TRUE)
  
  svg(file.path(opt$outdir, "goseq_pwf_fit.svg"),
      width = 800 / 96, height = 600 / 96)
  pwf <- nullp(gene_vec, bias.data = bias, plot.fit = TRUE)
  dev.off()
  
  go_res <- tryCatch(
    goseq(pwf, opt$genome, opt$gene_id_type, test.cats = c("GO:BP","GO:MF","GO:CC")),
    error = function(e) { message("goseq GO failed: ", e$message); NULL }
  )
  if (!is.null(go_res)) {
    go_res$padj <- p.adjust(go_res$over_represented_pvalue, method = "BH")
    write_tsv(go_res, file.path(opt$outdir, "goseq_GO_results.tsv"))
    
    hg_res <- goseq(pwf, opt$genome, opt$gene_id_type,
                    test.cats = "GO:BP", method = "Hypergeometric")
    hg_res$padj <- p.adjust(hg_res$over_represented_pvalue, method = "BH")
    
    merged <- inner_join(
      go_res |> filter(ontology == "BP") |>
        dplyr::select(category, term, goseq_p    = over_represented_pvalue,
               goseq_padj = padj),
      hg_res |> dplyr::select(category,       hyper_p    = over_represented_pvalue,
                       hyper_padj = padj),
      by = "category"
    )
    write_tsv(merged, file.path(opt$outdir, "goseq_vs_hypergeometric.tsv"))
    
    ggsave(file.path(opt$outdir, "goseq_vs_hypergeometric.svg"),
           ggplot(merged, aes(-log10(hyper_padj), -log10(goseq_padj))) +
             geom_point(alpha = 0.5) +
             geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
             labs(x = "-log10 adj.P (hypergeometric, biased)",
                  y = "-log10 adj.P (goseq, length-corrected)",
                  title = "Length-corrected vs uncorrected GO-BP"),
           width = 6.5, height = 6)
  }
  
} # end lengths_available (goseq)

message(">>> Done.")
