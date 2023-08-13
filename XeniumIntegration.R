library(Seurat)


wild.path <- "./data/Xenium_V1_FFPE_wildtype_13_4_months_outs"
disease.path <- "./data/Xenium_V1_FFPE_TgCRND8_17_9_months_outs"


wild.obj <- LoadXenium(wild.path, fov = "fov")
# remove cells with 0 counts
wild.obj <- subset(wild.obj, subset = nCount_Xenium > 0)

ImageDimPlot(wild.obj, fov = "fov", molecules = c("Gad1", "Sst", "Pvalb", "Gfap"), nmols = 20000)

ImageFeaturePlot(wild.obj, features = c("Cux2", "Rorb", "Bcl11b", "Foxp2"), max.cutoff = 'q90', size = 0.75, cols = c("white", "red"))

wild.obj <- SCTransform(wild.obj, assay = "Xenium")
wild.obj <- RunPCA(wild.obj, npcs = 30, features = rownames(wild.obj))
wild.obj <- RunUMAP(wild.obj, dims = 1:30)
wild.obj <- FindNeighbors(wild.obj, reduction = "pca", dims = 1:30)
wild.obj <- FindClusters(wild.obj, resolution = 0.3)

DimPlot(wild.obj)

FeaturePlot(wild.obj, features = c("Cux2", "Bcl11b", "Foxp2", "Gad1", "Sst", "Gfap"))

ImageDimPlot(wild.obj, cols = "polychrome", size = 0.75)


disease.obj <- LoadXenium(disease.path, fov = "fov")
# remove cells with 0 counts
disease.obj <- subset(disease.obj, subset = nCount_Xenium > 0)

ImageDimPlot(disease.obj, fov = "fov", molecules = c("Gad1", "Sst", "Pvalb", "Gfap"), nmols = 20000)

ImageFeaturePlot(disease.obj, features = c("Cux2", "Rorb", "Bcl11b", "Foxp2"), max.cutoff = 'q90', size = 0.75, cols = c("white", "red"))


disease.obj <- SCTransform(disease.obj, assay = "Xenium")
disease.obj <- RunPCA(disease.obj, npcs = 30, features = rownames(disease.obj))
disease.obj <- RunUMAP(disease.obj, dims = 1:30)
disease.obj <- FindNeighbors(disease.obj, reduction = "pca", dims = 1:30)
disease.obj <- FindClusters(disease.obj, resolution = 0.3)

DimPlot(disease.obj)

FeaturePlot(disease.obj, features = c("Cux2", "Bcl11b", "Foxp2", "Gad1", "Sst", "Gfap"))

ImageDimPlot(disease.obj, cols = "polychrome", size = 0.75)


#######################
#Integration
######################
wild.obj$disease <- "wild"
disease.obj$disease <- "TgCRND8"

combined.list <- list(wild.obj, disease.obj)

features <- SelectIntegrationFeatures(object.list = combined.list)

combined.list <- lapply(X = combined.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# select features that are repeatedly variable across datasets for integration
anchors <- FindIntegrationAnchors(object.list = combined.list, reduction = "rpca")

combined.obj <- IntegrateData(anchorset = anchors)


DefaultAssay(combined.obj) <- "integrated"

# Run the standard workflow for visualization and clustering
combined.obj <- ScaleData(combined.obj, verbose = FALSE)
combined.obj <- RunPCA(combined.obj, npcs = 30, verbose = FALSE)
combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:30)
combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:30)
combined.obj <- FindClusters(combined.obj, resolution = 0.5)

p1 <- DimPlot(combined.obj, reduction = "umap", group.by = "disease", raster=FALSE)
p2 <- DimPlot(combined.obj, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)
p1 + p2

DimPlot(combined.obj, reduction = "umap", split.by = "disease", raster=FALSE)

RNA<-combined.obj@assays$integrated
combined.obj@assays$RNA<-RNA


DefaultAssay(combined.obj) <- "RNA"
nk.markers <- FindConservedMarkers(combined.obj, ident.1 = 6, grouping.var = "disease", verbose = FALSE)
head(nk.markers)


FeaturePlot(combined.obj , features = c("Nxph3", "Plp1", "Slc17a7", "Apod", "Cpne6", "Nrn1"), min.cutoff = "q9", raster=FALSE)


mk <- FindMarkers(combined.obj, group.by = "disease", ident.1="wild", ident.2 = "TgCRND8",  verbose = FALSE)
head(mk, 15)

p1 <- ImageDimPlot(combined.obj, fov = "fov", molecules = c("Vip", "Nts", "Penk", "Pvalb"), nmols = 50000, alpha = 0.1)
p2 <- ImageDimPlot(combined.obj, fov = "fov.1", molecules = c("Vip", "Nts", "Penk", "Pvalb"), nmols = 50000, alpha = 0.1)
p1 + p2

> FeaturePlot(combined.obj , features = c("Vip", "Nts", "Penk", "Pvalb"), split.by = "disease", min.cutoff = "q9", raster=FALSE)