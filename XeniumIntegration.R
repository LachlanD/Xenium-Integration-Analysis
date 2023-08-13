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

RNA<-combined.obj@assays$integrated
combined.obj@assays$RNA<-RNA


DefaultAssay(combined.obj) <- "RNA"

# Run the standard workflow for visualization and clustering
combined.obj <- ScaleData(combined.obj, verbose = FALSE)
combined.obj <- RunPCA(combined.obj, npcs = 30, verbose = FALSE)
combined.obj <- RunUMAP(combined.obj, reduction = "pca", dims = 1:30)
combined.obj <- FindNeighbors(combined.obj, reduction = "pca", dims = 1:30)
combined.obj <- FindClusters(combined.obj, resolution = 0.5)

p1 <- DimPlot(combined.obj, reduction = "umap", group.by = "disease", raster=FALSE)
p2 <- DimPlot(combined.obj, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)
p1 + p2

DimPlot(combined.obj, reduction = "umap", split.by = "disease", raster=FALSE, label = TRUE)

p1 <- ImageDimPlot(combined.obj, fov = "fov.1",cols = "red", cells = WhichCells(combined.obj, idents = 2))
p2 <- ImageDimPlot(combined.obj, fov = "fov", ,cols = "red", cells = WhichCells(combined.obj, idents = 2))
p1 + p2


##############################################
#Find markers for annotating clusters
##############################################
nk.markers <- FindConservedMarkers(combined.obj, ident.1 = 2, grouping.var = "disease", verbose = FALSE)
head(nk.markers)

nk.markers <- FindConservedMarkers(combined.obj, ident.1 = 8, grouping.var = "disease", verbose = FALSE)
nk6 <- head(nk.markers, 6)


FeaturePlot(combined.obj , features = nk6, min.cutoff = "q9", raster=FALSE)


########################################################
#DGE between disease and wild within the same cluster
########################################################

combined.obj$celltype.disease <- paste(Idents(combined.obj), combined.obj$disease, sep = "_")
combined.obj$celltype <- Idents(combined.obj)
Idents(combined.obj) <- "celltype.disease"


#######################################
#Cluster 2 Analysis
#######################################
mk <- FindMarkers(combined.obj,  ident.1="2_wild", ident.2 = "2_TgCRND8",  verbose = FALSE)
top4 <- rownames(head(mk, 4))


p1 <- ImageDimPlot(combined.obj, fov = "fov.1", molecules = top4, nmols = 50000, cols = "white", alpha = 0.5, cells = WhichCells(combined.obj, idents = "2_TgCRND8"))
p2 <- ImageDimPlot(combined.obj, fov = "fov", molecules = top4, nmols = 50000, cols = "white", alpha = 0.5, cells = WhichCells(combined.obj, idents = "2_wild"))
p1 + p2

FeaturePlot(combined.obj , features = top4, split.by = "disease", min.cutoff = "q9", raster=FALSE)


#######################################
#Cluster 8 Analysis
#######################################
mk <- FindMarkers(combined.obj,  ident.1="8_wild", ident.2 = "8_TgCRND8",  verbose = FALSE)
top4 <- rownames(head(mk, 4))


p1 <- ImageDimPlot(combined.obj, fov = "fov.1", molecules = top4, nmols = 50000, cols = "white", alpha = 0.5, cells = WhichCells(combined.obj, idents = "8_TgCRND8"))
p2 <- ImageDimPlot(combined.obj, fov = "fov", molecules = top4, nmols = 50000, cols = "white", alpha = 0.5, cells = WhichCells(combined.obj, idents = "8_wild"))
p1 + p2

FeaturePlot(combined.obj , features = top4, split.by = "disease", min.cutoff = "q9", raster=FALSE)

