process CREATE {
    publishDir params.outdir, mode:'copy'

    input:
    path data_dir

    output:
    path 'pbmc3k.rds', emit: rds
    path 'VlnPlot.png', emit: png

    script:
    """#!/usr/bin/env Rscript

## Seurat
library(Seurat)
library(dplyr)

pbmc.data <- Read10X(data.dir = "${data_dir}")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
png(filename = "VlnPlot.png", res = 200, width = 1400, height = 600)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

saveRDS(pbmc, file = "pbmc3k.rds")

    """
} 