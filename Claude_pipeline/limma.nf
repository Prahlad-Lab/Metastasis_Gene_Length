process LIMMA_LENGTH_ANALYSIS {
    label 'process_medium'

    // params.limma_env_path is the full path to the pre-built conda environment
    // (supplied at runtime by submit_Metastasis_pipeline.sh via --limma_env_path).
    // Nextflow treats a string containing '/' as an existing environment path and
    // activates it directly without attempting to create a new environment.
    conda      params.limma_env_path
    container  null

    input:
    path expression
    path phenotype

    output:
    path 'results/*.tsv',              emit: tables
    path 'results/*.svg', optional: true, emit: plots
    path 'results/*.rds',              emit: rds

    script:
    def lengths_arg = params.gene_lengths ? "--gene_lengths ${file(params.gene_lengths)}" : ''
    """
    limma_length_analysis.R \\
        --expression    ${expression} \\
        --phenotype     ${phenotype} \\
        --data_type     ${params.data_type} \\
        --genome        ${params.genome_build} \\
        --gene_id_type  ${params.gene_id_type} \\
        --fdr           ${params.fdr_cutoff} \\
        --lfc           ${params.lfc_cutoff} \\
        ${lengths_arg} \\
        --outdir        results
    """
}
