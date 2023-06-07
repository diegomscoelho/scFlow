#!/usr/bin/env Rscript

## Seurat
library(Seurat)
library(dplyr)

pbmc <- readRDS("${rds}")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Top markers
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png(filename = "Top10genes.png", res = 200, width = 1600, height = 2000)
DoHeatmap(pbmc, features = top10\$gene) + NoLegend()
dev.off()

saveRDS(pbmc, file = "pbmc3k_final.rds")