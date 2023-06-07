process MARKERS {
    publishDir params.outdir, mode:'copy'

    input:
    path rds

    output:
    path 'pbmc3k_final.rds'
    path 'Top10genes.png'

    script:
    template 'find_markers.R'

} 