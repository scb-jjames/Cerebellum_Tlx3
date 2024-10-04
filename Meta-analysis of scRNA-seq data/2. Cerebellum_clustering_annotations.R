
#Find Tlx3 levels in different cell types of cerebellum at different dev stages
#Data is from : https://doi.org/10.1016/j.cub.2018.07.062
#Analysis by Rahul Jose
#

library(tidyverse)
library(gdata)
library(Seurat)
library(clusterProfiler)


# Clustering -------------------------------------------------------------------------

cerebx1 <- readRDS("cerebxs_QCed.rds")
cerebx2 <- NormalizeData(cerebx1)

cerebx2 <- FindVariableFeatures(cerebx2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(cerebx2)
cerebx2 <- ScaleData(cerebx2, features = all.genes)

saveRDS(object = cerebx2, file = "cerebx2_norm.scaled.rds")
cerebx2 <- RunPCA(cerebx2, features = VariableFeatures(object = cerebx2))

DimPlot(cerebx2, reduction = "pca")
ggsave(filename =  "PCA/PCA2_.jpeg",plot=last_plot(), path = "plots/", dpi = 320)
DimHeatmap(cerebx2, dims = 1:20, cells = 500, balanced = TRUE)
ggsave(filename =  "PCA/PCA3_.jpeg",plot=last_plot(), path = "plots/", dpi = 320)

cerebx2 <- ScoreJackStraw(object = cerebx2, dims = 1:20)
JackStrawPlot(cerebx2, dims = 1:20)
ggsave(filename =  "PCA/Jackstraw_.jpeg",plot=last_plot(), path = "plots/", dpi = 320)

ElbowPlot(cerebx2)
ggsave(filename =  "PCA/Elbow_.jpeg",plot=last_plot(), path = "plots/", dpi = 320)

cerebx2 <- FindNeighbors(cerebx2, dims = 1:20)
cerebx2 <- FindClusters(cerebx2, resolution = 0.5)

cerebx2 <- RunUMAP(object = cerebx2, dims = 1:20 )
tiff(filename = "UMAP_plot1.tiff",
     height = 5, width = 10,
     units = "in", res = 600)
DimPlot(cerebx2, reduction = "umap", label =F,
        split.by = "orig.ident")
dev.off()

saveRDS(object = cerebx2, file = "Clusters/cerebx2.rds")

# Cluster Annotation ------------------------------------------------------
all.markers <- FindAllMarkers(cerebx2, test.use = 'wilcox')
write.csv(x = all.markers , file = "all_markers_new.csv")

#1 # scMayoMap as a reference
library("scMayoMap")
gc()
scMayoMap.obj <- scMayoMap(data = all.markers, database=scMayoMapDatabase, tissue = 'brain')
pltx <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)
ggsave(filename = "clusters_newclustering_scMayoMap.jpeg", plot = pltx, dpi = 600)
scMayoMap2_markers <- scMayoMap.obj$markers
write.csv(x = scMayoMap2_markers, file = "scMayoMap_markers_new.csv")

#2 Markers from https://doi.org/10.1016/j.cub.2018.07.062 (except Tlx3) to be used

cereb.markers <- read.csv( file = "cerebellum_markers.csv")
cereb.markers <- dplyr::rename(cereb.markers, term = name)

head(rownames(cerebx2))
srt.genes <- rownames(cerebx2)

all.markers2 <- all.markers |> 
  filter(p_val_adj<=0.05 & avg_log2FC>0)

#enrich marker genes from literature on the clusters
enrich_clust <- function(avglfc_cutoff){
  x <- avglfc_cutoff
  
  for(i in unique(all.markers2$cluster)){
    
    tryCatch({
      
      paste("avg_log2FC>", x)
      paste("Cluster ", i)
      
      enrichment <- enricher(
        gene= filter(all.markers2, cluster == i)$gene,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 1,
        maxGSSize = 500,
        qvalueCutoff = 0.2,
        TERM2GENE=cereb.markers
      )
      en_res <- enrichment@result
      write.csv(en_res,
                file = paste("Clusters/plots/", "scb-atlas_lfc_cut-",x, "_cluster_", i,"_enrich_result.csv"))
      plt <- ggplot(en_res) +
        geom_point(aes(
          x = GeneRatio,
          y = Description,
          fill = p.adjust,
          color = p.adjust,
          size = Count
        ))
      ggsave(filename = paste("scb-atlas_lfc_cut-",x, "cluster_", i,"_enrich_result.jpeg"),
             plot = plt,
             path = "Clusters/plots/",
             dpi = 600
      )
      rm(enrichment, i, en_res, plt)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

enrich_clust(0)
enrich_clust(1)
enrich_clust(2)
enrich_clust(2.5)

# further distinguish GN and GNP clusters 
VlnPlot(cerebx2, features = c( "Mfap4", "Atoh1", "Neurod1", "Gabra6", "Pax6", "Lmx1a", "Math1", "Pde1c", "Pcsk9"), pt.size = 0)
ggsave("Clusters/plots/GN-GNP_subclusters.jpeg", plot=last_plot())


#Clusters -------------------------------------------------------------------------

# Rename clusters
new_names <- c("Granule neurons",              #0
               "Granule neuron progenitors",   #1
               "Progenitors",                  #2
               "Granule neuron progenitors",   #3
               "Purkinje cells",               #4
               "Granule neurons",              #5
               "Bergmann glia",                #6
               "Glial progenitor",             #7
               "Granule neuron progenitors",   #8
               "Glutamergic CN",               #9
               "GABAergic progenitors",        #10
               "Granule neuron progenitors",   #11
               "GABAergic interneurons",       #12
               "GABAergic interneurons",       #13
               "Midbrain cells",               #14
               "Endothelial cells",            #15
               "Roof plate cells",             #16
               "Microglia",                    #17
               "Meninges",                     #18
               "Meninges",                     #19
               "Oligodendrocytes",             #20
               "Erythrocytes",                 #21
               "Meninges",                     #22
               "Ciliated cells"                #23
)               

levels(cerebx2)
names(new_names) <- levels(cerebx2)
others <- c( "Midbrain cells", "Endothelial cells", "Roof plate cells",  "Microglia", "Meninges", "Erythrocytes", "Ciliated cells")
cerebcellstypes <- setdiff(new_names, others)

cerebx2$orig.clusters <- Idents(cerebx2)
cerebxcells5 <- Seurat::RenameIdents(object = cerebx2, new_names)

cerebxcells5 <- RunUMAP(object = cerebxcells5, dims = 1:20)

tiff(filename = "Clusters/UMAP_plot2.tiff",
     height = 5, width = 10,
     units = "in", res = 600)
DimPlot(cerebxcells5, reduction = "umap", label =T,
        split.by = "orig.ident")
dev.off()

saveRDS(object = cerebxcells5, file = "cerebxcells5_cluster-id.rds")
cerebxcells5 <- readRDS(file = "cerebxcells5_cluster-id.rds")

# -------------------------------------------------------------------------


