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

# ---------- CLI ----------
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
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load expression + phenotype
# =============================================================================
message(">>> Loading expression: ", opt$expression)
expr <- read_tsv(opt$expression, show_col_types = FALSE)
gene_col <- colnames(expr)[1]
mat <- expr |> column_to_rownames(var = gene_col) |> as.matrix()

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

# Paired design only makes sense if every patient appears in both groups
paired <- all(tapply(pheno$tissue_type, pheno$patient,
                     function(x) length(unique(x)) == 2))
message(">>> Design: ",
        if (paired) "paired (~patient + tissue_type)"
        else        "unpaired (~tissue_type)")

# =============================================================================
# 2. Branch on data type
# =============================================================================
use_voom <- opt$data_type %in% c("counts_raw", "counts_log2")

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
    dge  <- calcNormFactors(dge, method = "TMM")
}

# =============================================================================
# 3. Design and fit
# =============================================================================
design <- if (paired) model.matrix(~ patient + tissue_type, data = pheno)
          else        model.matrix(~ tissue_type,           data = pheno)

if (use_voom) {
    png(file.path(opt$outdir, "voom_mean_variance.png"),
        width = 800, height = 600, res = 120)
    v <- voomWithQualityWeights(dge, design, plot = TRUE)
    dev.off()
    fit <- lmFit(v, design)
} else {
    message(">>> Using limma-trend on log2-normalised expression")
    # Drop near-zero rows (log2-space noise floor)
    row_mean <- rowMeans(mat, na.rm = TRUE)
    keep <- row_mean > quantile(row_mean, 0.20, na.rm = TRUE)
    mat  <- mat[keep, ]
    message("    kept ", nrow(mat), " genes after noise-floor filter")
    v   <- mat   # for downstream uniformity
    fit <- lmFit(mat, design)
}

fit <- eBayes(fit, trend = !use_voom, robust = TRUE)

coef_name <- "tissue_typemetastasis"
stopifnot(coef_name %in% colnames(fit$coefficients))

res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none") |>
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
    len_vec <- getlength(res$gene_id, opt$genome, opt$gene_id_type)
    gene_len <- tibble(gene_id = res$gene_id, length_bp = len_vec)
}
res <- res |> left_join(gene_len, by = "gene_id")

write_tsv(res, file.path(opt$outdir, "limma_de_results.tsv"))
saveRDS(list(fit = fit, pheno = pheno, res = res, data_type = opt$data_type),
        file.path(opt$outdir, "limma_fit.rds"))

# =============================================================================
# 5. Length-stratified summaries + plots
# =============================================================================
res_len <- res |>
    filter(!is.na(length_bp), length_bp > 0) |>
    mutate(
        log10_len = log10(length_bp),
        decile    = ntile(length_bp, 10),
        sig       = adj.P.Val < opt$fdr,
        direction = case_when(sig & logFC >  opt$lfc ~ "up",
                              sig & logFC < -opt$lfc ~ "down",
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

theme_set(theme_bw(base_size = 11))

ggsave(file.path(opt$outdir, "length_vs_logfc.png"),
    ggplot(res_len, aes(log10_len, logFC)) +
        geom_hex(bins = 60) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_smooth(method = "loess", se = FALSE, colour = "red") +
        scale_fill_viridis_c(trans = "log10") +
        labs(x = "log10 gene length (bp)",
             y = "log2 FC (metastasis vs primary)",
             title = "Length vs fold change",
             subtitle = sprintf("Spearman rho = %.3f, p = %.2g",
                                sp$estimate, sp$p.value)),
    width = 7, height = 5, dpi = 150)

ggsave(file.path(opt$outdir, "length_decile_de_pct.png"),
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
    width = 8, height = 4.5, dpi = 150)

ggsave(file.path(opt$outdir, "length_density_by_direction.png"),
    ggplot(res_len, aes(log10_len, colour = direction, fill = direction)) +
        geom_density(alpha = 0.15) +
        scale_colour_manual(values = c(up = "#d73027", down = "#4575b4", ns = "grey50")) +
        scale_fill_manual  (values = c(up = "#d73027", down = "#4575b4", ns = "grey50")) +
        labs(x = "log10 gene length (bp)", y = "density",
             title = "Gene length distribution by DE direction"),
    width = 7, height = 4, dpi = 150)

ggsave(file.path(opt$outdir, "volcano.png"),
    ggplot(res_len, aes(logFC, -log10(adj.P.Val), colour = direction)) +
        geom_point(size = 0.6, alpha = 0.5) +
        geom_vline(xintercept = c(-opt$lfc, opt$lfc), linetype = 2) +
        geom_hline(yintercept = -log10(opt$fdr),       linetype = 2) +
        scale_colour_manual(values = c(up = "#d73027", down = "#4575b4", ns = "grey70")) +
        labs(title = "Volcano plot", x = "log2 FC", y = "-log10 adj.P"),
    width = 7, height = 5, dpi = 150)

# =============================================================================
# 6. goseq: length-bias-corrected GO/KEGG
# =============================================================================
message(">>> goseq enrichment")
gene_vec <- as.integer(res$adj.P.Val < opt$fdr & abs(res$logFC) > opt$lfc)
names(gene_vec) <- res$gene_id
gene_vec[is.na(gene_vec)] <- 0
bias <- setNames(res$length_bp, res$gene_id)
bias[is.na(bias)] <- median(bias, na.rm = TRUE)

png(file.path(opt$outdir, "goseq_pwf_fit.png"),
    width = 800, height = 600, res = 120)
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
                  select(category, term, goseq_p   = over_represented_pvalue,
                                          goseq_padj = padj),
        hg_res |> select(category,        hyper_p   = over_represented_pvalue,
                                          hyper_padj = padj),
        by = "category"
    )
    write_tsv(merged, file.path(opt$outdir, "goseq_vs_hypergeometric.tsv"))

    ggsave(file.path(opt$outdir, "goseq_vs_hypergeometric.png"),
        ggplot(merged, aes(-log10(hyper_padj), -log10(goseq_padj))) +
            geom_point(alpha = 0.5) +
            geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "red") +
            labs(x = "-log10 adj.P (hypergeometric, biased)",
                 y = "-log10 adj.P (goseq, length-corrected)",
                 title = "Length-corrected vs uncorrected GO-BP"),
        width = 6.5, height = 6, dpi = 150)
}

message(">>> Done.")
