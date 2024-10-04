# Single cell RNA-Seq of Tlx3 WT vs KO samples
# Budhaditya Basu
library(Seurat)
library(tidyverse)

# Create control seurat object
control.data <- Read10X("JJ_2_count/filtered_feature_bc_matrix/")
control <- CreateSeuratObject(counts = control.data,
                              project = "Control", 
                              min.cells = 3, 
                              min.features = 200)

# Check control meta data
control@meta.data

grep ("^mt-", rownames(control[["RNA"]]), value = T)

control[["Percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")

VlnPlot(control, features = c("nFeature_RNA", "nCount_RNA", 
                                    "Percent.mt"), ncol = 3)
control <- subset(control,
                  subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & Percent.mt < 7.5)

# Normalization of data
control <- NormalizeData(control)

# Identification of highly variable features (feature selection)
control <- FindVariableFeatures(control,
                                selection.method = "vst",
                                nfeatures = 2000)
#=============================================================================
# Create knock out seurat object
ko.data <- Read10X("KBH_count/filtered_feature_bc_matrix/")
ko <- CreateSeuratObject(counts = ko.data,
                         project = "KO", 
                         min.cells = 3, 
                         min.features = 200)

# Check KO meta data
ko@meta.data

ko[["Percent.mt"]] <- PercentageFeatureSet(ko, pattern = "^mt-")

VlnPlot(ko, features = c("nFeature_RNA", "nCount_RNA", 
                              "Percent.mt"), ncol = 3)
ko <- subset(ko,
             subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & Percent.mt < 7.5)

# Normalization of data
ko <- NormalizeData(ko)

# Identification of highly variable features (feature selection)

ko <- FindVariableFeatures(ko,
                           selection.method = "vst",
                           nfeatures = 2000)

#==============================================================================
# Data Integration
# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(control, ko))

tlx3.anchors <- FindIntegrationAnchors(object.list = list(control, ko),
                                       anchor.features = features)

tlx3.combined <- IntegrateData(anchorset = tlx3.anchors)




#==============================================================================
tiff(filename = "Data_QC.tiff", width = 10, height = 5, units = "in", res = 600)
VlnPlot(tlx3.combined, features = c("nFeature_RNA", "nCount_RNA", 
                                    "Percent.mt"), ncol = 3, 
        split.by = "orig.ident")
dev.off()

# scaling for PCA
all.genes <- rownames(tlx3.combined)
tlx3.combined <- ScaleData(tlx3.combined, features = all.genes)

# Perform linear dimensional reduction
tlx3.combined <- RunPCA(tlx3.combined, 
                        features = VariableFeatures(object = tlx3.combined))

#Elbow plot 

tiff(filename = "elbow_plot_tlx3.tiff", width = 7, height = 5, 
     units = "in",
     res = 600)
ElbowPlot(tlx3.combined)
dev.off()


#============================================================

#Cluster the cells finalized 
tlx3.combined <- FindNeighbors(tlx3.combined, dims = 1:20)

tlx3.combined <- FindClusters(tlx3.combined, resolution = 0.4)

tlx3.combined <- RunUMAP(tlx3.combined, dims = 1:20)

tlx3.combined@meta.data

##Plot the UMAP plot split by orig.ident
tiff(filename = "UMAP_plot_PN5_split.tiff",
     height = 5, width = 10,
     units = "in", res = 600)
DimPlot(tlx3.combined, reduction = "umap", label = TRUE,
        split.by = "orig.ident")+
  NoLegend()
dev.off()

# Grouped UMAP plot grouped by orig.ident
tiff(filename = "UMAP_plot_grouped.tiff", 
     height = 6, width = 8,
     units = "in", res = 600)
DimPlot(tlx3.combined, group.by = "orig.ident")+labs(title = "")
dev.off()

#Save as RDS file
saveRDS(tlx3.combined, file = "sc_PN5_tlx3.combined_integrated.rds")

# Check Cre expression and Atoh1 expression
DefaultAssay(tlx3.combined)
DefaultAssay(tlx3.combined) <- "RNA"

FeaturePlot(tlx3.combined, features = c("Cre", "Atoh1"),
            split.by = "orig.ident")

