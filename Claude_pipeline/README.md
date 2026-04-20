# metastasis-length-pipeline (Xena edition)

Nextflow DSL2 pipeline that pulls primary vs metastasis RNA-seq from
**UCSC Xena** via `UCSCXenaTools`, runs `limma` DE, and does explicit
gene-length bias analysis (Mandelboum 2019 / Modur 2018 framework).

## Quickstart

```bash
# TCGA SKCM (melanoma) — primary vs metastatic. ~100 vs ~360 samples.
nextflow run main.nf -profile conda

# Toil recompute (TCGA+TARGET+GTEx unified), filtered to melanoma
nextflow run main.nf -profile conda \
    --preset pancan_toil \
    --cancer_types 'Skin Cutaneous Melanoma'

# Custom: any Xena dataset + phenotype
nextflow run main.nf -profile conda \
    --preset custom \
    --xena_host   toilHub \
    --xena_dataset   'TcgaTargetGtex_RSEM_Hugo_norm_count' \
    --xena_phenotype 'TcgaTargetGTEX_phenotype.txt' \
    --data_type log2_norm
```

## Presets

| preset         | source                      | data type     | notes                                       |
|----------------|-----------------------------|---------------|---------------------------------------------|
| `skcm`         | TCGA SKCM HiSeqV2           | `log2_norm`   | ~100 primary + ~360 metastatic (default)    |
| `pancan_toil`  | Toil TCGA+TARGET+GTEx       | `log2_norm`   | Uniformly reprocessed, use `--cancer_types` |
| `custom`       | user-specified host+dataset | user          | Specify all three `--xena_*` flags          |

### Data type → which limma variant

- `counts_raw` — rare on Xena; uses `voomWithQualityWeights`.
- `counts_log2` — e.g. Toil's `TcgaTargetGtex_RSEM_gene_expected_count`
  stored as `log2(x+1)`. Reverse-transformed, rounded, then `voom`.
- `log2_norm` — e.g. TCGA HiSeqV2 `log2(norm_count+1)`. Uses
  `limma` with `eBayes(trend=TRUE)` directly; `voom` is not appropriate
  because the data are already library-normalised in log space.

## Sample filtering

Samples are partitioned on TCGA barcode positions 14–15:

| code  | meaning                    | default group     |
|-------|----------------------------|-------------------|
| `01`  | Primary Solid Tumor        | primary           |
| `02`  | Recurrent Solid Tumor      | (excluded)        |
| `06`  | Metastatic                 | metastasis        |
| `07`  | Additional Metastatic      | metastasis        |
| `11`  | Solid Tissue Normal        | (excluded)        |

Override with `--primary_codes` / `--metastasis_codes`.

Patient ID is derived from TCGA barcode positions 1–12, so a paired
design is detected automatically when a patient has both a primary
and a metastatic sample (rare in TCGA, more common in MET500).

## Output

```
results/
├── xena_download/
│   ├── expression.tsv.gz
│   ├── phenotype.tsv
│   └── xena_download_log.txt
└── limma_length_analysis/
    ├── limma_de_results.tsv            # topTable + gene length
    ├── length_decile_summary.tsv
    ├── length_logfc_correlation.tsv
    ├── goseq_GO_results.tsv
    ├── goseq_vs_hypergeometric.tsv
    ├── voom_mean_variance.png          # (counts modes only)
    ├── length_vs_logfc.png             # main length-bias plot
    ├── length_decile_de_pct.png
    ├── length_density_by_direction.png
    ├── volcano.png
    ├── goseq_pwf_fit.png
    ├── goseq_vs_hypergeometric.png
    └── limma_fit.rds
```

## Notes / caveats

- **Network**: `XENA_DOWNLOAD` needs outbound HTTPS to the relevant
  Xena hub. On HPCs without internet, pre-download and point the
  pipeline at local files (I can add a flag for that).
- **Gene ID types**: SKCM HiSeqV2 and Toil-Hugo use HGNC symbols, so
  use `--gene_id_type geneSymbol`. Ensembl-keyed datasets need
  `--gene_id_type ensGene`. Version suffixes (`.N`) are stripped
  automatically.
- **Length source**: by default, `goseq::getlength()` provides lengths
  from the UCSC knownGene/ensGene tracks matching `--genome_build`.
  For custom annotation, pass `--gene_lengths your_lengths.tsv`.
- **Sample size for SKCM preset**: while primary-vs-metastatic
  comparison is well-powered (~100 vs 360), most primary samples
  are excised lesions that may already harbour invasive programs.
  Interpret as "metastatic deposit vs invasive primary" rather than
  "early primary vs late metastasis".

## References

- Wang S. *et al* (2019) UCSCXenaTools: an R package for accessing
  UCSC Xena data. *JOSS* 4(40):1627.
- Goldman M.J. *et al* (2020) Visualizing and interpreting cancer
  genomics data via the Xena platform. *Nat Biotechnol* 38:675.
- Vivian J. *et al* (2017) Toil enables reproducible, open source,
  big biomedical data analyses. *Nat Biotechnol* 35:314. (Toil recompute)
- Robinson D.R. *et al* (2017) Integrative clinical genomics of
  metastatic cancer. *Nature* 548:297. (MET500)
- Modur V. *et al* (2018) Defective transcription elongation in a
  subset of cancers confers immunotherapy resistance. *Nat Commun* 9:4410.
- Mandelboum S. *et al* (2019) Recurrent functional misinterpretation
  of RNA-seq data caused by sample-specific gene length bias.
  *PLOS Biol* 17:e3000481.
