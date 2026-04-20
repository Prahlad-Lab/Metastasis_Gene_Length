process LIMMA_LENGTH_ANALYSIS {
    label 'process_medium'

    conda 'bioconda::bioconductor-limma=3.58.1 bioconda::bioconductor-edger=4.0.16 bioconda::bioconductor-goseq=1.54.0 conda-forge::r-optparse conda-forge::r-tidyverse conda-forge::r-ggrepel'
    container '../meta_xena_limma.sif'

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
