#!/usr/bin/env nextflow
/*
 * metastasis-length-pipeline (Xena edition)
 *
 * Pulls RNA-seq expression from UCSC Xena (default: TCGA SKCM or the
 * Toil recompute hub) and runs limma-based DE + gene-length bias analysis
 * for primary vs metastatic samples.
 */

nextflow.enable.dsl = 2

include { XENA_DOWNLOAD         } from './xena_download.nf'
include { LIMMA_LENGTH_ANALYSIS } from './limma.nf'

workflow {

    XENA_DOWNLOAD()

    LIMMA_LENGTH_ANALYSIS(
        XENA_DOWNLOAD.out.expression,
        XENA_DOWNLOAD.out.phenotype
    )
}

workflow.onComplete {
    log.info "Done. Results: ${params.outdir}"
}
