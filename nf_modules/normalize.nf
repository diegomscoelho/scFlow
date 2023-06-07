process NORMALIZE {
    publishDir params.outdir, mode:'copy'

    input:
    path rds

    output:
    path 'pbmc3k_normalized.rds', emit: rds
    path 'DimensionPlot.png'

    script:
    """#!/usr/bin/env Rscript

## Seurat
library(Seurat)
library(dplyr)

pbmc <- readRDS("${rds}")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Create Dimension Plots
pca <- DimPlot(pbmc, reduction = "pca")
umap <- DimPlot(pbmc, reduction = "umap")

png(filename = "DimensionPlot.png", res = 200, width = 1600, height = 1000)
pca + umap
dev.off()

saveRDS(pbmc, file = "pbmc3k_normalized.rds")

    """
} 