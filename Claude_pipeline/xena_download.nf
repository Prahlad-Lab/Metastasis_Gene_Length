process XENA_DOWNLOAD {
    label 'process_low'

    conda 'conda-forge::r-ucscxenatools conda-forge::r-tidyverse conda-forge::r-optparse'
    container '/vscratch/grp-vprahlad/Metastasis_Gene_Length/meta_xena.sif'

    output:
    path 'expression.tsv.gz', emit: expression
    path 'phenotype.tsv',     emit: phenotype
    path 'xena_download_log.txt'

    script:
    def types  = params.cancer_types ? "--cancer_types '${params.cancer_types}'" : ''
    """
    xena_download.R \\
        --host       ${params.xena_host} \\
        --dataset    ${params.xena_dataset} \\
        --phenotype  ${params.xena_phenotype} \\
        --primary_codes    ${params.primary_codes} \\
        --metastasis_codes ${params.metastasis_codes} \\
        ${types} \\
        --outdir . 2>&1 | tee xena_download_log.txt
    """
}
