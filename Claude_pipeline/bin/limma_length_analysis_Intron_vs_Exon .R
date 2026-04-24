#!/usr/bin/env Rscript

# =============================================================================
# limma DE + gene-length bias analysis, matrix-input (Xena) version.
# Structural Update: Added Intron vs. Exon (CDS) Length Burden Analysis
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(limma)
  library(edgeR)
  library(goseq)
  library(ggplot2)
  library(ggpubr)    # Added for stat_cor() in scatter plots
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(biomaRt)
})

# ---------------------------------------------------------------------------
# Parameter Setup (Handles both RStudio Interactive and CLI Execution)
# ---------------------------------------------------------------------------
if (interactive()) {
  # --- RSTUDIO INTERACTIVE MODE ---
  opt <- list(
    expression   = "/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/results/xena_download/expression.tsv.gz", 
    phenotype    = "/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/results/xena_download/phenotype.tsv",    
    data_type    = "log2_norm",                         
    sample_col   = "sample",                            
    group_col    = "sample_type",                       
    baseline     = "Primary Tumor",                     
    contrast     = "Metastatic",                        
    genome       = "hsapiens_gene_ensembl",             
    gene_id_type = "hgnc_symbol",                       
    outdir       = "/vscratch/grp-vprahlad/Metastasis_Gene_Length/Claude_pipeline/results_intron_exon"                
  )
} else {
  # --- CLI MODE ---
  option_list <- list(
    make_option(c("-e", "--expression"), type="character", help="Path to expression matrix TSV/CSV."),
    make_option(c("-p", "--phenotype"), type="character", help="Path to phenotype TSV/CSV."),
    make_option(c("-d", "--data_type"), type="character", default="counts_raw", help="Data type: 'counts_raw', 'counts_log2', or 'log2_norm'."),
    make_option(c("-s", "--sample_col"), type="character", default="sample", help="Column name in phenotype for sample IDs."),
    make_option(c("-g", "--group_col"), type="character", default="sample_type", help="Column name in phenotype for group labels."),
    make_option(c("-b", "--baseline"), type="character", help="Baseline group name (e.g., 'Primary')."),
    make_option(c("-c", "--contrast"), type="character", help="Contrast group name (e.g., 'Metastasis')."),
    make_option(c("--genome"), type="character", default="hsapiens_gene_ensembl", help="biomaRt dataset (default: hsapiens_gene_ensembl)."),
    make_option(c("--gene_id_type"), type="character", default="hgnc_symbol", help="Gene ID type (e.g., hgnc_symbol, ensembl_gene_id)."),
    make_option(c("-o", "--outdir"), type="character", default="results", help="Output directory.")
  )
  opt <- parse_args(OptionParser(option_list=option_list))
}

# Ensure output directory exists
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

message("======================================================")
message("limma DE & Gene Length Analysis (Structural Edition)")
message("======================================================")
for (n in names(opt)) {
  message(sprintf("  %-15s : %s", n, opt[[n]]))
}
message("------------------------------------------------------\\n")

# 1. Load Expression Data
message("Loading expression data...")
expr <- read_tsv(opt$expression, show_col_types = FALSE)

gene_col <- colnames(expr)[1]
message(sprintf("Assuming column 1 ('%s') contains gene IDs.", gene_col))

expr_mat <- expr |>
  column_to_rownames(var = gene_col) |>
  as.matrix()

# Handle potential duplicate gene names
if (any(duplicated(expr[[gene_col]]))) {
  message("WARNING: Duplicate gene IDs found. Aggregating by summing...")
  expr_mat <- expr |>
    group_by(!!sym(gene_col)) |>
    summarise(across(everything(), sum, na.rm = TRUE)) |>
    column_to_rownames(var = gene_col) |>
    as.matrix()
}

# 2. Load Phenotype Data
message("Loading phenotype data...")
pheno <- read_tsv(opt$phenotype, show_col_types = FALSE)

if (!opt$sample_col %in% colnames(pheno)) stop(sprintf("Sample column '%s' not found in phenotype data.", opt$sample_col))
if (!opt$group_col %in% colnames(pheno)) stop(sprintf("Group column '%s' not found in phenotype data.", opt$group_col))

pheno <- pheno |>
  filter(!!sym(opt$group_col) %in% c(opt$baseline, opt$contrast)) |>
  mutate(group = factor(!!sym(opt$group_col), levels = c(opt$baseline, opt$contrast)))

common_samples <- intersect(colnames(expr_mat), pheno[[opt$sample_col]])
message(sprintf("Found %d samples in common between expression and phenotype.", length(common_samples)))

expr_mat <- expr_mat[, common_samples, drop=FALSE]
pheno <- pheno |> filter(!!sym(opt$sample_col) %in% common_samples) |> arrange(match(!!sym(opt$sample_col), common_samples))

# 3. Transform Data Based on Type
message(sprintf("Processing data as '%s'...", opt$data_type))
if (opt$data_type == "counts_raw") {
  counts_mat <- round(expr_mat)
} else if (opt$data_type == "counts_log2") {
  counts_mat <- round(2^expr_mat - 1)
  counts_mat[counts_mat < 0] <- 0
} else if (opt$data_type == "log2_norm") {
  norm_mat <- expr_mat 
} else {
  stop("Unknown data_type. Use 'counts_raw', 'counts_log2', or 'log2_norm'.")
}

# 4. Fetch Gene Lengths (Transcript, CDS, and Intron)
message("Fetching Transcript and CDS lengths from Ensembl...")

gene_data <- NULL

# Helper to run the query safely
run_bm_query <- function(mart_obj) {
  tryCatch({
    getBM(
      attributes = c(opt$gene_id_type, 'transcript_length', 'cds_length', 'transcript_biotype'),
      mart = mart_obj
    )
  }, error = function(e) { return(NULL) })
}

# Attempt 1: US East Mirror
message("  -> Attempting US East Mirror...")
try({
  mart <- useEnsembl(biomart = "ensembl", dataset = opt$genome, mirror = "useast")
  gene_data <- run_bm_query(mart)
}, silent = TRUE)

# Attempt 2: US West Mirror
if (is.null(gene_data)) {
  message("  -> Attempting US West Mirror...")
  try({
    mart <- useEnsembl(biomart = "ensembl", dataset = opt$genome, mirror = "uswest")
    gene_data <- run_bm_query(mart)
  }, silent = TRUE)
}

# Attempt 3: Default Live Server
if (is.null(gene_data)) {
  message("  -> Attempting Default Live Ensembl Server...")
  try({
    mart <- useMart("ensembl", dataset = opt$genome)
    gene_data <- run_bm_query(mart)
  }, silent = TRUE)
}

# Attempt 4: The "Failsafe" Archive (Ensembl Release 109 - Feb 2023)
if (is.null(gene_data)) {
  message("  -> Live servers down. Attempting Stable Ensembl Archive (Feb 2023)...")
  try({
    mart <- useMart(host="https://feb2023.archive.ensembl.org", 
                    biomart="ENSEMBL_MART_ENSEMBL", 
                    dataset=opt$genome)
    gene_data <- run_bm_query(mart)
  }, silent = TRUE)
}

# Final Check
if (is.null(gene_data)) {
  stop("CRITICAL ERROR: All Ensembl servers (Mirrors, Main, and Archives) are currently unreachable. Ensembl is likely experiencing a severe global outage. Please try running the script again in an hour.")
}

message("Successfully retrieved length data! Processing Intron Burden...")
gene_lengths <- gene_data |>
  filter(transcript_biotype == "protein_coding") |>
  filter(!!sym(opt$gene_id_type) != "") |>
  mutate(cds_length = ifelse(is.na(cds_length), 0, cds_length)) |>
  group_by(!!sym(opt$gene_id_type)) |>
  arrange(desc(transcript_length)) |>
  slice(1) |>
  ungroup() |>
  mutate(intron_length = transcript_length - cds_length) |>
  dplyr::select(
    !!sym(opt$gene_id_type), 
    transcript_length, 
    cds_length, 
    intron_length
  ) |>
  distinct()

# 5. Differential Expression via limma
message("Running limma...")
design <- model.matrix(~ group, data = pheno)

if (opt$data_type %in% c("counts_raw", "counts_log2")) {
  dge <- DGEList(counts = counts_mat)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  v <- voomWithQualityWeights(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
} else {
  fit <- lmFit(norm_mat, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE) 
}

if (opt$data_type %in% c("counts_raw", "counts_log2")) {
  fit <- eBayes(fit)
}

res <- topTable(fit, coef = 2, number = Inf, sort.by = "none") |>
  rownames_to_column("ID")

write_tsv(res, file.path(opt$outdir, "limma_results_full.tsv"))

# 6. Merge DE Results with Gene Lengths
message("Merging DE results with gene lengths...")
res_length <- res |>
  inner_join(gene_lengths, by = c("ID" = opt$gene_id_type)) |>
  na.omit()

# 7. Structural Analysis: CDS vs Intron Burden
message("======================================================")
message("Performing Structural Intron vs. CDS Analysis...")
message("======================================================")

# Filter for significant genes to calculate the structural bias models
sig_res <- res_length |> filter(adj.P.Val < 0.05)

if (nrow(sig_res) > 10) {
  sig_res <- sig_res |>
    mutate(
      log10_CDS = log10(cds_length + 1),
      log10_Intron = log10(intron_length + 1),
      Significance = ifelse(logFC > 0, "Up", "Down")
    )
  
  # Multiple Linear Regression
  message("Running Multiple Linear Regression (logFC ~ log10_CDS + log10_Intron)...")
  lm_model <- lm(logFC ~ log10_CDS + log10_Intron, data = sig_res)
  model_summary <- summary(lm_model)
  
  sink(file.path(opt$outdir, "Structural_Regression_Stats.txt"))
  cat("=== Multiple Regression Analysis: Metastasis vs Primary ===\n")
  cat("Formula: logFC ~ log10(CDS Length) + log10(Intron Length)\n\n")
  print(model_summary)
  sink()
  
  pval_cds <- formatC(model_summary$coefficients["log10_CDS", "Pr(>|t|)"], format = "e", digits = 2)
  pval_intron <- formatC(model_summary$coefficients["log10_Intron", "Pr(>|t|)"], format = "e", digits = 2)
  
  # Plot 1: CDS Length
  p_cds <- ggplot(sig_res, aes(x = log10_CDS, y = logFC)) +
    geom_point(aes(color = Significance), alpha = 0.4, size = 1) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red")) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", formula = 'y ~ x') +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(title = "FC vs CDS Length",
         subtitle = paste("Regression p-val:", pval_cds),
         x = "Log10(CDS Length [bp])", y = "Log2 Fold Change") +
    theme_minimal() + theme(legend.position = "none")
  
  ggsave(file.path(opt$outdir, "Corr_CDS_Length.pdf"), p_cds, width = 6, height = 5)
  
  # Plot 2: Intron Length
  p_intron <- ggplot(sig_res, aes(x = log10_Intron, y = logFC)) +
    geom_point(aes(color = Significance), alpha = 0.4, size = 1) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red")) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", formula = 'y ~ x') +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    labs(title = "FC vs Intron Length",
         subtitle = paste("Regression p-val:", pval_intron),
         x = "Log10(Total Intron Length [bp])", y = "Log2 Fold Change") +
    theme_minimal() + theme(legend.position = "none")
  
  ggsave(file.path(opt$outdir, "Corr_Intron_Length.pdf"), p_intron, width = 6, height = 5)
  
} else {
  message("Not enough significant DEGs to run structural regression modeling.")
}


# 8. Run goseq (Length-Bias Corrected Enrichment)
message("Running goseq (Length-bias corrected Enrichment)...")
genes_vector <- as.integer(res_length$adj.P.Val < 0.05)
names(genes_vector) <- res_length$ID

if (sum(genes_vector) > 0) {
  pwf <- nullp(genes_vector, bias.data = res_length$transcript_length, plot.fit = FALSE)
  
  go_res <- tryCatch(
    suppressMessages(goseq(pwf, opt$genome, opt$gene_id_type, test.cats = c("GO:BP"))),
    error = function(e) { message("goseq GO failed: ", e$message); NULL }
  )
  
  if (!is.null(go_res)) {
    go_res$padj <- p.adjust(go_res$over_represented_pvalue, method = "BH")
    write_tsv(go_res, file.path(opt$outdir, "goseq_GO_results.tsv"))
    message("goseq completed successfully.")
  }
} else {
  message("No significant genes found (adj.P.Val < 0.05). Skipping goseq.")
}

message("======================================================")
message("Pipeline Finished. Outputs saved in: ", normalizePath(opt$outdir))
message("======================================================")