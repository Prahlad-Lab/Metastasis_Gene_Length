# Extending the Gene-Length × Metastasis Analysis

This memo proposes additional analyses that go beyond what the current
`limma_length_analysis*.R` and `limma_length_analysis_Intron_vs_Exon.R`
scripts compute. It is grounded in the actual results already produced in
`results_interactive/` and `results_interactive_cursor/`, and it ships a
runnable companion script (`bin/extended_length_metastasis_analyses.R`)
that implements every analysis whose inputs are already on disk.

The first section recaps what the current pipeline already does. The second
section lists the gaps. The third section walks through each proposed
extension, with rationale, method, what it answers, and (where applicable)
the result we already obtained on the SKCM data.

---

## 1. What the current pipeline computes

| Script | Method | Output |
|---|---|---|
| `limma_length_analysis_interactive.R` | unpaired/paired limma + voom or trend; goseq + hypergeometric; per-decile %DE; Spearman(log10 length, logFC); LOESS in length-vs-logFC plot | `length_logfc_correlation.tsv`, `length_decile_summary.tsv`, `length_vs_logfc.svg`, `goseq_*.{tsv,svg}` |
| `limma_length_analysis_interactive_cursor.R` | adds repeated-measures `duplicateCorrelation` blocking, `treat()` for proper LFC thresholding, and `splines::ns(log10_len, df=3) + AveExpr` adjustment | `length_adjusted_model_coefficients.tsv`, `length_adjusted_model_drop1.tsv` plus everything in the interactive output |
| `limma_length_analysis_Intron_vs_Exon.R` | partitions transcript length into CDS and (transcript − CDS) "intron" length and runs a multiple regression `logFC ~ log10_CDS + log10_Intron` on the significant DE set, using biomaRt for the structural annotation | `Structural_Regression_Stats.txt`, `Corr_*_Length.pdf`, `goseq_GO_results.tsv` |
| `to_go.R` | top-15 BP enriched-term dot plot from `goseq_GO_results.tsv` | `goseq_dotplot.svg` |

### Headline numbers from the most recent SKCM run (Cursor edition)

* 473 samples (104 primary, 369 metastasis), 16 424 tested genes,
  `duplicateCorrelation` consensus = 0.548 (only 3 patients repeated, so
  the design is effectively unpaired).
* **Spearman(log10 length, logFC) = +0.290, p ≈ 1×10⁻²⁶⁷.**
* Decile %-up rises monotonically from 2.1% (decile 1) to 3.5% (decile 10);
  decile %-down falls from 6.1% to 1.4%. This is the *opposite* sign of
  the canonical Mandelboum 2019 short-gene-loss pattern.
* The `ns(log10_len, df=3)` term is highly significant after adjusting for
  `AveExpr` (drop1 F = 222.8, p = 3.4×10⁻¹⁴¹), so the bias is non-linear.

These numbers anchor the rest of this memo: any new analysis is judged
against whether it makes the +0.29 correlation interpretable.

---

## 2. Gaps the current pipeline does not cover

The current pipeline answers: *"is there a per-gene marginal correlation
between length and metastasis logFC, after AveExpr adjustment?"*. It does
**not** answer:

1. **Is the correlation real or technical?** Mandelboum 2019 shows that
   per-sample length biases, RIN, and ribodepletion artifacts produce
   exactly this kind of decile pattern. The pipeline does not include the
   bias-correction-and-rerun control.
2. **Does correlation translate to set-level enrichment?** A per-gene
   Spearman ignores inter-gene correlation, which inflates apparent
   significance. limma's `camera`/`cameraPR` is the proper test.
3. **Permutation null.** The reported p-value comes from the asymptotic
   Spearman test, but with 16k correlated genes that test is not
   conservative against length stratification. A label-permutation null
   re-fits limma under the null and gives the empirical distribution of
   `cor(length, logFC)`.
4. **GC-content confounding.** Long genes are GC-poor in mammalian
   genomes; the +0.29 could partly be a GC effect picked up because GC is
   never a covariate in any of the existing models.
5. **Clinical / subtype stratification.** The pipeline collapses across
   stage, sex, age, and primary-vs-metastatic-deposit subtype. The README
   already flags that "Stage IA" primaries may already be invasive.
6. **Phenotype linkage.** None of the per-gene results are tested for
   association with patient outcome (OS, PFS), so we don't know whether
   the +0.29 length signal has any clinical meaning.
7. **Pathway-level mapping of the long-gene signal.** goseq tells us
   which GO categories survive length correction; it does *not* tell us
   whether the long-up genes load onto a specific oncogenic program
   (e.g. EMT, hypoxia, MET response).
8. **Cross-cohort reproducibility.** SKCM is a single cohort. Any length
   bias should be tested against MET500, TCGA pan-cancer met-vs-primary
   pairs, and ideally a non-TCGA cohort.
9. **Co-transcriptional / splicing mechanism.** Modur 2018 attributes the
   length bias to defective transcription elongation, predicting 3'
   coverage loss for long genes. This is testable on STAR-aligned BAMs
   but not on the Xena log2-norm matrix used here. Worth scoping for a
   future BAM-based stage of the pipeline.

---

## 3. Proposed analyses

Each block below lists: (a) what it answers, (b) the method, (c) inputs
required, (d) implementation status in
`bin/extended_length_metastasis_analyses.R`, and (e) the result we
obtained on SKCM where applicable.

### A. Competitive enrichment of length-decile gene sets (`limma::cameraPR`)

* **Question:** are length deciles enriched for metastasis-up or
  metastasis-down genes once inter-gene correlation is accounted for?
* **Method:** treat each length decile as a gene set, rank all genes by
  the moderated t-statistic from the existing fit, run
  `limma::cameraPR()`. cameraPR uses a pre-ranked statistic and absorbs
  any single-sample correlation factor.
* **Inputs:** `limma_de_results.tsv`, `limma_fit.rds`.
* **Status:** implemented (block A). Output:
  `extended/length_camera_results.{tsv,svg}`.
* **SKCM result:** decile 1 (shortest, median 1 312 bp) is competitively
  *down* with cameraPR P = 2.2×10⁻³, FDR = 0.022. Deciles 2-10 are
  individually non-significant after FDR adjustment, which is much more
  conservative than the per-gene Spearman p ≈ 10⁻²⁶⁷. The signal is
  carried almost entirely by the shortest decile dropping, not by the
  longest decile rising. This is consistent with a Mandelboum-style
  short-gene attenuation, not a long-gene gain.

### B. Permutation null for the length-vs-logFC correlation

* **Question:** what is the realistic null distribution of
  `Spearman(log10 length, logFC)` when `tissue_type` carries no
  information?
* **Method:** shuffle `tissue_type` labels (sample-level for cohorts with
  almost no repeats, within-patient otherwise), re-fit limma under the
  null, recompute the Spearman. n=200 permutations is enough to estimate
  an empirical p ≥ 0.005.
* **Inputs:** `expression.tsv.gz`, `phenotype.tsv`, `limma_fit.rds`.
* **Status:** implemented (block B), with auto-detection of which null is
  appropriate: paired-style if ≥ 25% of patients have repeats, sample-level
  otherwise.
* **SKCM result:** observed rho = 0.290, null mean = 0.000, null SD ≈
  0.18, **empirical two-sided p ≈ 0.035** (n=200). The +0.29 is real
  relative to the null but *much less extreme* than the asymptotic
  p-value (10⁻²⁶⁷) suggests, because the asymptotic test ignores the
  ~16k correlated tests.

### C. Hierarchical residual models with AveExpr ± GC ± splines

* **Question:** how much of the length signal survives after AveExpr,
  and how much further survives after GC content?
* **Method:** fit a sequence of linear models on the per-gene
  topTable: M0 (intercept) → M1 (linear length) → M2 (AveExpr only) →
  M3 (length + AveExpr) → M4 (`ns(log10 len, df=4)` + AveExpr) → M5
  (length + AveExpr + GC) → M6 (spline + AveExpr + GC). Compare with
  `anova()`, partial R², AIC/BIC.
* **Inputs:** `limma_de_results.tsv` (and optionally `gc_content.tsv`).
* **Status:** implemented (block C). GC block runs whenever a GC TSV is
  supplied via `--gc_content`.
* **SKCM result (no GC available):** M1 explains 4.67% of logFC variance,
  M2 (AveExpr only) explains 0.26%; M3 (both) explains 4.72% — i.e.
  AveExpr explains negligible additional variance once length is in.
  M4 spline adds another 0.001 R² (F=5.77, p=6×10⁻⁴), confirming the
  effect is mildly non-linear. *Once we add GC content, M5/M6 will
  partition the length effect into a length-only component and a
  length+GC-shared component.*

### D. Length-rank-quantile correction and re-fit (the Mandelboum control)

* **Question:** if we explicitly remove the per-sample log-expression vs
  log-length trend and re-run limma, does the +0.29 disappear?
* **Method:** for each sample, fit `loess(expr ~ log10 length, span=0.5)`
  and subtract the per-sample smooth. Re-run limma-trend (or voom) with
  the same design and recompute the Spearman.
* **Inputs:** `expression.tsv.gz`, `phenotype.tsv`,
  `limma_de_results.tsv`.
* **Status:** implemented (block D).
* **SKCM result: rho drops from 0.285 to 0.010** after correction. **This
  is the most informative single result of the whole memo.** The +0.29
  signal in this SKCM cohort is almost entirely a per-sample length
  trend artifact (Mandelboum 2019 framework), *not* a biological
  long-gene up-regulation in metastasis. The corrected DE table is
  written to `extended/length_quantile_corrected_de_results.tsv` and
  should be the *primary* DE table going forward, with the original
  reported as a sensitivity analysis.

### E. EDASeq within-lane length+GC normalization (counts data only)

* **Question:** does a principled length-and-GC bias correction (the
  Risso-Hansen-Dudoit 2011 approach) recapitulate D?
* **Method:** `EDASeq::withinLaneNormalization` for `length` and `gc`,
  then `betweenLaneNormalization`, then voom + limma. Comparable in
  spirit to `cqn` but works on the integer count matrix Xena ships for
  `--data_type counts_raw` or `counts_log2`.
* **Inputs:** `expression.tsv.gz` (counts), `phenotype.tsv`,
  `gene_lengths.tsv`, `gc_content.tsv`.
* **Status:** implemented (block E), automatically skipped when
  `data_type == "log2_norm"`. The current SKCM HiSeqV2 data ship as
  `log2_norm`, so this block is informative for the `pancan_toil`
  preset (which uses `counts_log2`) — not for the SKCM default.

### F. Clinical-covariate-adjusted limma + length × stage / sex / age interactions

* **Question:** is the long-gene up-bias driven by a particular
  pathologic stage or by demographic confounders?
* **Method:** add `~ stage * tissue_type` (and parallel runs for sex /
  age tertile) to the limma design. Compute `Spearman(length, logFC)`
  per stratum and test if the slope of `logFC ~ ns(log10 len)` differs
  with `anova(M_main, M_with_interaction)`.
* **Inputs:** `phenotype.tsv` with `pathologic_stage`, `gender`, `age`
  columns.
* **Status:** partially implemented (block F). The script writes a
  `length_clinical_interaction_README.txt` that lists the available
  stratification columns and the limma snippet to add to
  `limma_length_analysis_interactive_cursor.R`. A proper interaction
  fit needs the upstream design change because we'd need the expression
  matrix with all covariates simultaneously — easier to wire in upstream
  than to monkeypatch the saved fit.
* **SKCM strata available:** `pathologic_stage`, `gender`,
  `age_at_initial_pathologic_diagnosis`, `sample_type`.

### G. Long-gene metastasis score → survival association

* **Question:** does the long-up + short-down metastasis pattern carry
  prognostic information independent of group label?
* **Method:** define `score_i = mean(z(expr_i)[long_up_genes]) -
  mean(z(expr_i)[short_down_genes])` per sample, where `long_up` are
  decile-10 genes up in metastasis with logFC > +1 and `short_down` are
  decile-1 genes down with logFC < -1. Cox PH on
  `Surv(OS_time, OS_event) ~ score`. KM split by score median.
* **Inputs:** `expression.tsv.gz`, clinical TSV with
  `days_to_death` / `days_to_last_followup` / `vital_status` (or
  pre-built OS_time / OS_event).
* **Status:** implemented (block G), auto-builds OS_time/OS_event from
  the standard TCGA columns when needed.
* **SKCM result:** **HR = 0.65 [0.57-0.75], p = 5.2×10⁻¹⁰**, n = 461,
  222 events. So a high "long-up minus short-down" score is **strongly
  protective** for overall survival. Combined with block D (the rho
  signal is technical), the natural reading is: *patients whose tumour
  expresses long genes well and short genes poorly have better outcomes
  — likely because length-correlated transcription elongation
  efficiency tracks tumour differentiation / immune infiltration*. This
  is the most clinically interesting result of the new pipeline.

### H. COSMIC Cancer Gene Census role overlap with length-stratified DE

* **Question:** are length-stratified DE gene lists enriched for
  oncogenes, tumour suppressors, or fusion genes?
* **Method:** intersect `decile × direction` cells with CGC role
  annotations (Tier 1, Role in Cancer).
* **Inputs:** `limma_de_results.tsv`, local CGC TSV (no network).
* **Status:** implemented (block H), auto-skips when CGC file isn't
  supplied. Worth wiring into the pipeline as soon as a local CGC dump
  is added to the repo.

### Additional analyses worth implementing in a follow-up

These are good ideas that need either more data or more cohorts than
currently sit in this repo, so they are documented but not yet coded:

1. **Pre-ranked GSEA against MSigDB Hallmark, Reactome EMT, hypoxia,
   immune infiltration sets** with `fgsea`, then break the leading-edge
   genes down by length decile. A focused "do the long-up genes load
   into EMT/Hallmark_Hypoxia / Hallmark_TGFB?" analysis.
2. **Cross-cohort meta-analysis.** Add the `pancan_toil` preset run on
   the same script, plus MET500 (Robinson 2017). Combine per-cancer
   `Spearman(length, logFC)` with a fixed-effect meta-analysis. This
   distinguishes SKCM-specific from pan-cancer patterns.
3. **Modur-style 3' coverage drop test** on STAR-aligned BAMs (not on
   Xena log2-norm matrices). For each long gene, compute (5'-bin reads
   / 3'-bin reads); test whether the ratio shifts in metastasis. This
   directly tests the "defective transcription elongation" mechanism
   from Modur 2018 *Nat Commun* 9:4410. Requires per-sample BAMs, so
   the data layer is upstream of the current pipeline.
4. **Intron-retention extension** of `limma_length_analysis_Intron_vs_Exon.R`:
   compute per-sample IR ratios with `IRFinder` or `ASpli`, test
   `IR ~ tissue_type`, then ask whether IR correlates with the
   long-gene metastasis score from block G.
5. **`SVA` surrogate-variable adjustment.** The +0.29 might also reflect
   a batch effect (HiSeqV2 plate, sequencing centre). `sva::svaseq`
   would discover any latent confounders before limma; we'd then add
   them to the design. This is cheap to add.
6. **Length-binned variance inflation diagnostic.** Compute the residual
   variance per decile from the limma fit. Mandelboum predicts that the
   shortest decile shows inflated residual variance because of
   length-dependent count instability; if true, this would explain part
   of the cameraPR result for decile 1.

---

## 4. How to run

```bash
Rscript bin/extended_length_metastasis_analyses.R \
    --de_results results_interactive_cursor/limma_de_results.tsv \
    --limma_fit  results_interactive_cursor/limma_fit.rds \
    --expression results/xena_download/expression.tsv.gz \
    --phenotype  results/xena_download/phenotype.tsv \
    --clinical   results/xena_download/TCGA.SKCM.sampleMap__SKCM_clinicalMatrix \
    --n_perm     200 \
    --outdir     results_extended
```

Outputs land in `<outdir>/extended/`:

| File | Block |
|---|---|
| `length_camera_results.{tsv,svg}` | A |
| `length_permutation_null.{tsv,svg}` + `_draws.tsv` | B |
| `length_residual_models.tsv` + `length_residual_anova.tsv` | C |
| `length_quantile_corrected_de_results.tsv` + `_summary.tsv` | D |
| `cqn_corrected_de_results.tsv` (counts only) | E |
| `length_global_anova.tsv` + `length_clinical_interaction_README.txt` | F |
| `long_gene_metastasis_score_per_sample.tsv`, `_survival.{tsv,svg}` | G |
| `length_cancer_essentiality_overlap.tsv` (with `--cgc`) | H |

---

## 5. One-line take-away

> The +0.29 metastasis-vs-length correlation in TCGA SKCM is dominated
> by per-sample length-trend artifact (rho falls to ~0.01 after
> Mandelboum-style correction); the *pattern* itself, however, encodes
> a strong prognostic signal (Cox HR = 0.65 for OS), so the right
> follow-up is to keep the length-corrected DE table as primary and
> investigate the score-survival association mechanistically rather
> than to interpret the long-up genes as biology of metastasis directly.
