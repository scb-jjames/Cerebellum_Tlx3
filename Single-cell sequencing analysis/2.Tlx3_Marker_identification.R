# Single cell RNA-Seq of Tlx3 WT vs KO samples
# Budhaditya Basu
library(Seurat)
library(tidyverse)
# Load RDS file
tlx3.combined <- readRDS("sc_PN5_tlx3.combined_integrated.rds")

DefaultAssay(tlx3.combined) <- "RNA"

# Check the orig.ident names
unique(tlx3.combined@meta.data$orig.ident)

# Calculate the proportion of each cell type comparing the conditions
# How many cells are in each cluster
table(Idents(tlx3.combined))

# How many cells are in each replicate?
table(tlx3.combined$orig.ident)
# Subset the seurat object



# Identify marker genes for all clusters
markers_seur <- Seurat::FindAllMarkers(tlx3.combined, 
                                       only.pos = TRUE,
                                       test.use = "wilcox",
                                       logfc.threshold = 0.25)
head(markers_seur)
# Saving marker genes to file
write.csv(markers_seur, 
          file = "combined_all_markers.csv", 
          quote = FALSE, 
          row.names = FALSE)



# Retrieve the top 5 marker genes per cluster
top5 <- markers_seur %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_seur), value = TRUE)),
                   n = 5)

# Create the dot plot
tiff(filename = "dot_plot_top5_marker_Tlx3.tiff", height = 5, 
     width = 10,
     units = "in", res = 600)
Seurat::DotPlot(tlx3.combined, features = unique(top5$gene)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                     size = 8, hjust = 1)) +
  Seurat::NoLegend()
dev.off()


DefaultAssay(tlx3.combined) <- "integrated"
#plot heatmap
tiff(filename = "top5_heatmap_all_cluster.tiff", height = 6, width = 12,
     units = "in", res = 600)
DoHeatmap(tlx3.combined, features = top5$gene) + NoLegend()
dev.off()

DefaultAssay(tlx3.combined) <- "RNA"

# Marker Identification


# Granule Cell Progenitor- Mfap4, Atoh1, Mki67, Mycn, Pax6, Tlx3 (Cluster 0 and Cluster 1)
#==============================================================================
VlnPlot(tlx3.combined, features = c("Atoh1", "Mfap4", "Cre"),
        split.by = "orig.ident")
FeaturePlot(tlx3.combined, features = c("Atoh1", "Mfap4", "Cre"),
            split.by = "orig.ident")

# Granule Cells --Neurod1, Cntn2, Reln, Pax6, Tlx3 (Cluster 2)
#=============================================================================
VlnPlot(tlx3.combined, features = c("Neurod1", "Cntn2", "Reln"),
        split.by = "orig.ident")
FeaturePlot(tlx3.combined, features = c("Neurod1", "Cntn2", "Reln"),
            split.by = "orig.ident")

# GABArgic Interneuron- Pax2 (Cluster 3)
#=========================================================================
VlnPlot(tlx3.combined, features = c("Pax2"),
        split.by = "orig.ident")
FeaturePlot(tlx3.combined, features = c("Pax2"),
            split.by = "orig.ident")

# Purkinje Cell - Calb1, Car8 (Cluster 4)
#===============================================================================
tiff(filename = "Violin_plot_Car8&Calb1.tiff",
     height = 6, width = 12,
     units = "in", res = 600)
VlnPlot(tlx3.combined, features = c("Car8", "Calb1"),
        split.by = "orig.ident")
dev.off()

# Bergmann glia-- Hopx and Gdf10 (Cluster 5)
#===============================================================================
##Astrocytes marker Slc1a3
FeaturePlot(tlx3.combined, features = c("Slc1a3"),
            split.by = "orig.ident")
VlnPlot(tlx3.combined, features = "Slc1a3", split.by = "orig.ident")

VlnPlot(tlx3.combined, features = c("Hopx", "Gdf10"),
        split.by = "orig.ident")
FeaturePlot(tlx3.combined, features = c("Hopx",  "Gdf10"),
            split.by = "orig.ident")

#Oligodendrocytes -- Sox10, Olig1 (Cluster 6)
#===============================================================================
VlnPlot(tlx3.combined, features = c("Sox10", "Olig1"),
        split.by = "orig.ident")
FeaturePlot(tlx3.combined, features = c("Sox10", "Olig1"),
            split.by = "orig.ident")

#Endothelial Precursor Cell -- Lgals1, Col3a1 (cluster 7)
#===============================================================================

VlnPlot(tlx3.combined, features = c("Lgals1", "Col3a1"),
        split.by = "orig.ident")
FeaturePlot(tlx3.combined, features = c("Lgals1", "Col3a1"),
            split.by = "orig.ident")

#Erythrocytes --Hba-a2, Alas2, Hbb-bt, Fech (Cluster 8)
#===============================================================================
VlnPlot(tlx3.combined, features = c("Hba-a2", "Alas2","Hbb-bt", "Fech"),
        split.by = "orig.ident")
FeaturePlot(tlx3.combined, features = "Hba-a2", split.by = "orig.ident")

#Endothelial cell markers: Cldn5, Sparc (Cluster 9)
#===============================================================================
VlnPlot(tlx3.combined, features = c("Cldn5", "Kdr"),
        split.by = "orig.ident")
FeaturePlot(tlx3.combined, features = c("Cldn5", "Kdr"), split.by = "orig.ident")

#================================================================================

#Atoh1_Tlx3_Cre expression
tiff(filename = "Atoh1_Tlx3_Cre_expression.tiff",
     height = 12, width = 12,
     units = "in", res = 600)
FeaturePlot(tlx3.combined, features = c("Atoh1", "Cre", "Tlx3"),
            pt.size = 2,
            split.by = "orig.ident")
dev.off()

#Violin Plot
tiff(filename = "Atoh1_Tlx3_Cre_expression_violin.tiff",
     height = 3, width = 10,
     units = "in", res = 600)
VlnPlot(tlx3.combined, features = c("Atoh1", "Tlx3", "Cre"),
        split.by = "orig.ident")
dev.off()

#Violin Plot
tiff(filename = "Atoh1_expression_violin.tiff",
     height = 3, width = 10,
     units = "in", res = 600)
VlnPlot(tlx3.combined, features = c("Atoh1"),
        split.by = "orig.ident")
dev.off()

#Violin Plot
tiff(filename = "Cre_expression_violin.tiff",
     height = 3, width = 10,
     units = "in", res = 600)
VlnPlot(tlx3.combined, features = c("Cre"),
        split.by = "orig.ident")
dev.off()

#Violin Plot
tiff(filename = "Tlx3_expression_violin.tiff",
     height = 3, width = 10,
     units = "in", res = 600)
VlnPlot(tlx3.combined, features = c("Tlx3"),
        split.by = "orig.ident")
dev.off()
#=================================================================================
cells.plot0 <- colnames(tlx3.combined[, 
                                     (tlx3.combined$seurat_clusters %in% c("0"))])

library(SCpubr)

tiff(filename = "Tlx3_expression_in_0_cluster.tiff",
     height = 6, width = 10,
     units = "in",
     res = 600)
SCpubr::do_FeaturePlot(sample = tlx3.combined, 
                       features = "Tlx3",
                       split.by = "orig.ident",
                       cells.highlight = cells.plot0,
                       pt.size = 0.2)
dev.off()
tlx3.combined@meta.data


###Annotate clusters based on marker identification
new.cluster.ids <- c("Granule Neuron Prog", 
                     "Granule Neuron Prog", 
                     "Granule Neuron", 
                     "GABAergic Interneuron",
                     "Purkinje Cells",
                     "Bergmann Glia",
                     "Oligodendrocytes",
                     "Endothelial Precursor Cell",
                     "Erythrocytes", "Endothelial cells")
names(new.cluster.ids) <- levels(tlx3.combined)
tlx3.combined <- RenameIdents(tlx3.combined, new.cluster.ids)

tlx3.combined$CellType <- Idents(tlx3.combined)
#Save Annotated Seurat Object
saveRDS(tlx3.combined, file = "sc_PN5_tlx3.combined_annotated.rds")




