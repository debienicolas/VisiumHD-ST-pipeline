
input_path <- "results/human_breast_cancer_final/bin2cell/bin2cell_output.h5ad"


data <- schard::h5ad2seurat_spatial(input_path)

# normalize data
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) # make some filtering based on QC metrics visualizations, see Seurat tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
data <- ScaleData(data, features = rownames(data))
data <- RunPCA(data, features = VariableFeatures(object = data))

# Check number of PC components (we selected 10 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(data)

# cluster and visualize
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.8)
data <- RunUMAP(data, dims = 1:10)
DimPlot(data, reduction = "umap")