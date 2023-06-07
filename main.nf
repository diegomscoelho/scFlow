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

include { CREATE    } from "$baseDir/nf_modules/create.nf"
include { NORMALIZE } from "$baseDir/nf_modules/normalize.nf"
include { MARKERS   } from "$baseDir/nf_modules/find_markers.nf"

workflow {
    create_ch =     CREATE(params.data_dir)
    normalize_ch =  NORMALIZE(create_ch.rds)
    markers_ch =    MARKERS(normalize_ch.rds)
 }

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir\n" : "Oops .. something went wrong" )
}
