#!/usr/bin/env Rscript

# =============================================================================
# Download expression + phenotype from UCSC Xena and emit:
#   expression.tsv.gz  (genes x samples, wide)
#   phenotype.tsv      (sample, patient, tissue_type, cancer_type, ...)
# =============================================================================

suppressPackageStartupMessages({
    library(UCSCXenaTools)
    library(optparse)
    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
    library(tibble)
})

# ---------- CLI ----------
option_list <- list(
    make_option("--host",             type = "character"),
    make_option("--dataset",          type = "character"),
    make_option("--phenotype",        type = "character"),
    make_option("--primary_codes",    type = "character", default = "01"),
    make_option("--metastasis_codes", type = "character", default = "06,07"),
    make_option("--cancer_types",     type = "character", default = NA),
    make_option("--outdir",           type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

primary_codes <- strsplit(opt$primary_codes,    ",")[[1]]
met_codes     <- strsplit(opt$metastasis_codes, ",")[[1]]

# =============================================================================
# 1. Find the Xena dataset record
# =============================================================================
message(">>> Looking up Xena dataset")
message("    host:      ", opt$host)
message("    dataset:   ", opt$dataset)
message("    phenotype: ", opt$phenotype)

data(XenaData, package = "UCSCXenaTools")

# Validate dataset exists in XenaData metadata
host_map <- c(tcgaHub = "https://tcga.xenahubs.net",
              toilHub = "https://toil.xenahubs.net",
              gdcHub  = "https://gdc.xenahubs.net",
              pancanAtlasHub = "https://pancanatlas.xenahubs.net",
              publicHub = "https://ucscpublic.xenahubs.net")

host_url <- host_map[opt$host]
if (is.na(host_url)) stop("Unknown --host: ", opt$host,
                          ". See names(host_map) in this script.")

# =============================================================================
# 2. Download expression matrix
# =============================================================================
message(">>> Downloading expression matrix")
expr_query <- XenaGenerate(subset = XenaHostNames == opt$host) |>
    XenaFilter(filterDatasets = paste0("^", opt$dataset, "$"))

if (nrow(expr_query@datasets) == 0) {
    stop("Dataset not found in Xena registry: ", opt$dataset)
}

expr_dl <- XenaQuery(expr_query) |>
    XenaDownload(destdir = opt$outdir, trans_slash = TRUE, force = FALSE)
expr <- XenaPrepare(expr_dl)   # data.frame: sample x genes OR genes x samples

# Xena returns data with a 'sample' column (=gene ids) and one column per sample.
# Rename for clarity.
if (colnames(expr)[1] != "sample") {
    colnames(expr)[1] <- "sample"
}
expr <- expr |> rename(gene = sample)
message("    expression dim: ", nrow(expr), " genes x ", ncol(expr) - 1, " samples")

# =============================================================================
# 3. Download phenotype
# =============================================================================
message(">>> Downloading phenotype")
pheno_query <- XenaGenerate(subset = XenaHostNames == opt$host) |>
    XenaFilter(filterDatasets = paste0("^", opt$phenotype, "$"))

if (nrow(pheno_query@datasets) == 0) {
    stop("Phenotype dataset not found: ", opt$phenotype)
}

pheno_dl <- XenaQuery(pheno_query) |>
    XenaDownload(destdir = opt$outdir, trans_slash = TRUE, force = FALSE)
pheno <- XenaPrepare(pheno_dl) |> as_tibble()

# Normalise the sample-id column name.
sample_col <- intersect(c("sampleID", "sample", "SAMPLE_ID"), colnames(pheno))[1]
if (is.na(sample_col)) stop("No sample id column found in phenotype")
pheno <- pheno |> rename(sample = !!sample_col)

# =============================================================================
# 4. Derive tissue_type from TCGA sample barcode (positions 14-15)
#    TCGA barcode: TCGA-XX-XXXX-##A-... where ## is the sample-type code.
# =============================================================================
extract_sample_type <- function(x) {
    # works for TCGA barcodes and TARGET (similar structure)
    code <- str_sub(x, 14, 15)
    code[!str_detect(x, "^TCGA-")] <- NA
    code
}

pheno <- pheno |>
    mutate(
        sample_type_code = extract_sample_type(sample),
        tissue_type = case_when(
            sample_type_code %in% primary_codes ~ "primary",
            sample_type_code %in% met_codes     ~ "metastasis",
            TRUE                                ~ NA_character_
        ),
        patient = str_sub(sample, 1, 12)      # TCGA patient barcode
    )

# Optional cancer-type filter (useful for pancan_toil preset)
if (!is.na(opt$cancer_types)) {
    keep_types <- strsplit(opt$cancer_types, ",")[[1]] |> str_trim()
    type_col <- intersect(
        c("primary_disease", "disease_code", "_primary_disease",
          "cancer type abbreviation", "detailed_category"),
        colnames(pheno)
    )[1]
    if (!is.na(type_col)) {
        pheno <- pheno |> filter(.data[[type_col]] %in% keep_types)
        message("    filtered by ", type_col, " in {",
                paste(keep_types, collapse = ", "), "}: ",
                nrow(pheno), " samples")
    } else {
        warning("Could not find a cancer-type column to filter on")
    }
}

# Keep only samples we have a tissue_type for AND that exist in the expression matrix
expr_samples <- setdiff(colnames(expr), "gene")
pheno <- pheno |>
    filter(!is.na(tissue_type), sample %in% expr_samples) |>
    select(sample, patient, tissue_type, sample_type_code, any_of(c(
        "primary_disease", "disease_code", "_primary_disease",
        "cancer type abbreviation", "gender", "age_at_initial_pathologic_diagnosis",
        "pathologic_stage", "sample_type"
    ))) |>
    distinct(sample, .keep_all = TRUE)

message(">>> Sample counts by tissue_type:")
print(table(pheno$tissue_type))

if (sum(pheno$tissue_type == "primary")    < 3 ||
    sum(pheno$tissue_type == "metastasis") < 3) {
    stop("Need at least 3 samples per group. Check cohort + sample-type codes.")
}

# Subset expression matrix to the kept samples
expr_out <- expr |> select(gene, all_of(pheno$sample))

# =============================================================================
# 5. Write outputs
# =============================================================================
write_tsv(expr_out, file.path(opt$outdir, "expression.tsv.gz"))
write_tsv(pheno,    file.path(opt$outdir, "phenotype.tsv"))
message(">>> Wrote expression.tsv.gz and phenotype.tsv to ", opt$outdir)
