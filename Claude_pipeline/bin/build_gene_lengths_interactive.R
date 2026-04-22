#!/usr/bin/env Rscript
# =============================================================================
# build_gene_lengths.R
#
# Pre-build a gene_id <TAB> length_bp TSV that can be supplied to the pipeline
# via --gene_lengths, avoiding runtime goseq/biomaRt network failures.
#
# Usage (run once, outside Nextflow):
#
#   Rscript build_gene_lengths.R \
#       --gene_id_type geneSymbol \
#       --outfile gene_lengths_hg38_symbol.tsv
#
#   Rscript build_gene_lengths.R \
#       --gene_id_type ensGene \
#       --outfile gene_lengths_hg38_ensg.tsv
#
# Then run the pipeline with:
#   nextflow run main.nf --gene_lengths /path/to/gene_lengths_hg38_symbol.tsv ...
#
# Strategy:
#   Uses TxDb.Hsapiens.UCSC.hg38.knownGene + org.Hs.eg.db (Bioconductor,
#   fully local — no network required) to compute the longest transcript length
#   per gene. Falls back to biomaRt if TxDb packages are not installed.
# =============================================================================

# ---------------------------------------------------------------------------
# Parameter Setup (Handles both RStudio Interactive and CLI Execution)
# ---------------------------------------------------------------------------
if (interactive()) {
  # --- RSTUDIO INTERACTIVE MODE ---
  # Manually set your parameters here when running line-by-line
  opt <- list(
    gene_id_type = "geneSymbol", # Options: "geneSymbol" or "ensGene"
    outfile      = "gene_lengths_interactive_test.tsv"
  )
  message(">>> Running interactively in RStudio. Using hardcoded parameters.")
} else {
  # --- COMMAND LINE MODE ---
  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option("--gene_id_type", type = "character", default = "geneSymbol",
                help = "Gene ID type: geneSymbol (default) or ensGene"),
    make_option("--outfile", type = "character", default = "gene_lengths.tsv",
                help = "Output TSV path (gene_id TAB length_bp)")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
}

message(">>> Building gene lengths for hg38 / ", opt$gene_id_type)

# ---------------------------------------------------------------------------
# Try local TxDb approach first (fastest, no network)
# ---------------------------------------------------------------------------
local_ok <- tryCatch({
  if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
    stop("TxDb.Hsapiens.UCSC.hg38.knownGene not installed")
  if (!requireNamespace("GenomicFeatures", quietly = TRUE))
    stop("GenomicFeatures not installed")
  
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  tx_by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
  # Compute transcript widths and take max per Entrez gene
  tx_widths <- GenomicFeatures::transcriptLengths(txdb, with.cds_len = FALSE)
  len_entrez <- tapply(tx_widths$tx_len, tx_widths$gene_id, max, na.rm = TRUE)
  
  if (opt$gene_id_type == "geneSymbol") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("org.Hs.eg.db not installed")
    eg2sym <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = names(len_entrez),
      columns = "SYMBOL",
      keytype = "ENTREZID"
    )
    eg2sym <- eg2sym[!is.na(eg2sym$SYMBOL) & !duplicated(eg2sym$ENTREZID), ]
    matched <- len_entrez[eg2sym$ENTREZID]
    out <- data.frame(gene_id   = eg2sym$SYMBOL,
                      length_bp = as.integer(matched),
                      stringsAsFactors = FALSE)
    out <- out[!is.na(out$length_bp), ]
    # Collapse duplicate symbols: keep max length
    out <- out[order(-out$length_bp), ]
    out <- out[!duplicated(out$gene_id), ]
  } else {
    # ensGene: map Entrez -> Ensembl
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("org.Hs.eg.db not installed")
    eg2ens <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = names(len_entrez),
      columns = "ENSEMBL",
      keytype = "ENTREZID"
    )
    eg2ens <- eg2ens[!is.na(eg2ens$ENSEMBL) & !duplicated(eg2ens$ENTREZID), ]
    matched <- len_entrez[eg2ens$ENTREZID]
    out <- data.frame(gene_id   = eg2ens$ENSEMBL,
                      length_bp = as.integer(matched),
                      stringsAsFactors = FALSE)
    out <- out[!is.na(out$length_bp), ]
    out <- out[order(-out$length_bp), ]
    out <- out[!duplicated(out$gene_id), ]
  }
  
  write.table(out, opt$outfile,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(">>> Wrote ", nrow(out), " entries to ", opt$outfile,
          " (via TxDb.Hsapiens.UCSC.hg38.knownGene)")
  TRUE
}, error = function(e) {
  message("    TxDb approach failed: ", e$message)
  message(">>> Falling back to biomaRt")
  FALSE
})

# Proceed to fallbacks only if local approach failed (avoids calling quit() in RStudio)
if (!local_ok) {
  
  # ---------------------------------------------------------------------------
  # biomaRt fallback (tries multiple Ensembl mirrors)
  # ---------------------------------------------------------------------------
  biomart_ok <- tryCatch({
    if (!requireNamespace("biomaRt", quietly = TRUE))
      stop("biomaRt not installed")
    
    # useEnsembl() is preferred over useMart() and supports mirror selection
    mart <- NULL
    mirrors <- c("useast", "uswest", "asia", "www")
    for (mirror in mirrors) {
      mart <- tryCatch(
        biomaRt::useEnsembl("ensembl",
                            dataset = "hsapiens_gene_ensembl",
                            mirror  = mirror),
        error = function(e) {
          message("    biomaRt mirror '", mirror, "' failed: ", e$message)
          NULL
        }
      )
      if (!is.null(mart)) {
        message("    biomaRt connected via mirror: ", mirror)
        break
      }
    }
    if (is.null(mart))
      stop("Could not connect to any Ensembl mirror")
    
    attrs <- if (opt$gene_id_type == "geneSymbol") {
      c("hgnc_symbol", "transcript_length")
    } else {
      c("ensembl_gene_id", "transcript_length")
    }
    
    bm <- biomaRt::getBM(attributes = attrs, mart = mart)
    bm <- bm[bm[[1]] != "" & !is.na(bm[[1]]), ]
    
    agg <- tapply(bm[[2]], bm[[1]], max, na.rm = TRUE)
    out <- data.frame(gene_id   = names(agg),
                      length_bp = as.integer(agg),
                      stringsAsFactors = FALSE)
    
    write.table(out, opt$outfile,
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    message(">>> Wrote ", nrow(out), " entries to ", opt$outfile, " (via biomaRt)")
    TRUE
  }, error = function(e) {
    message("    biomaRt approach failed: ", e$message)
    message(">>> Falling back to geneLenDataBase")
    FALSE
  })
  
  # Proceed to final fallback only if biomaRt failed
  if (!biomart_ok) {
    # ---------------------------------------------------------------------------
    # geneLenDataBase fallback (local, no network required)
    # ---------------------------------------------------------------------------
    tryCatch({
      if (!requireNamespace("geneLenDataBase", quietly = TRUE))
        stop("geneLenDataBase not installed")
      
      data_name <- paste0("hg38.", opt$gene_id_type, ".LENGTH")
      tryCatch(
        utils::data(list = data_name, package = "geneLenDataBase"),
        error = function(e) stop("Dataset '", data_name, "' not found in geneLenDataBase")
      )
      len_data <- get(data_name)
      
      out <- data.frame(gene_id   = names(len_data),
                        length_bp = as.integer(len_data),
                        stringsAsFactors = FALSE)
      out <- out[!is.na(out$length_bp) & out$length_bp > 0, ]
      
      write.table(out, opt$outfile,
                  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      message(">>> Wrote ", nrow(out), " entries to ", opt$outfile,
              " (via geneLenDataBase)")
    }, error = function(e) {
      stop("All gene length retrieval methods failed: ", e$message)
    })
  }
}
