/*
 * pipeline input parameters
 */

// params.barcodes = "$projectDir/data/barcodes.tsv"
// params.matrix = "$projectDir/data/matrix.mtx"
// params.features = "$projectDir/data/genes.tsv"
params.data_dir = "$projectDir/data/"
params.outdir = "results"

// log.info """\
//     S E U R A T - N F   P I P E L I N E
//     ===================================
//     Barcodes        : ${params.barcodes}
//     Matrix          : ${params.matrix}
//     Features        : ${params.features}
//     """
//     .stripIndent()

include { CREATE_SEURAT } from "$baseDir/nf_modules/create.nf"

workflow {
    create_ch = CREATE_SEURAT(params.data_dir))
 }

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir\n" : "Oops .. something went wrong" )
}
