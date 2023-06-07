process MARKERS {
    publishDir params.outdir, mode:'copy'

    input:
    path rds

    output:
    path 'pbmc3k_final.rds'

    script:
    """#!/usr/bin/env Rscript

## Seurat
library(Seurat)
library(dplyr)

pbmc <- readRDS("${rds}")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top markers
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(pbmc, features = top10\$gene) + NoLegend()

saveRDS(pbmc, file = "pbmc3k_final.rds")

    """
} 