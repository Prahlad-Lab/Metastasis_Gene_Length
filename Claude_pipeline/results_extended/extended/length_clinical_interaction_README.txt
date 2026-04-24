Phenotype columns available for stratified re-fits:
  pathologic_stage
  gender
  age_at_initial_pathologic_diagnosis
  sample_type

To run a true tissue_type * <stratum> interaction limma fit, re-run
the upstream limma_length_analysis_interactive_cursor.R after adding the
interaction term to its design matrix, e.g.:
    design <- model.matrix(~ pathologic_stage * tissue_type, data = pheno)
and pass --outdir results_extended_<stratum>.
