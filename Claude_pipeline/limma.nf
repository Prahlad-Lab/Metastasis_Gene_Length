process LIMMA_LENGTH_ANALYSIS {
    label 'process_medium'

    conda      'test_limma_env'
    container  null

    input:
    path expression
    path phenotype

    output:
    path 'results/*.tsv', emit: tables
    path 'results/*.svg', emit: plots
    path 'results/*.rds', emit: rds

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
