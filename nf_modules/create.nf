process CREATE {
    publishDir params.outdir, mode:'copy'

    input:
    path data_dir

    output:
    path 'pbmc3k.rds', emit: rds
    path 'VlnPlot.png', emit: png

    script:
    template 'create.R'
} 