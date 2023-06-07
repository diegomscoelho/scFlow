process NORMALIZE {
    publishDir params.outdir, mode:'copy'

    input:
    path rds

    output:
    path 'pbmc3k_normalized.rds', emit: rds
    path 'DimensionPlot.png'

    script:
    template "normalize.R"
    
} 