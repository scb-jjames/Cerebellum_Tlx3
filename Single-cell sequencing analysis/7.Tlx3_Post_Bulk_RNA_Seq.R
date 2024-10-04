# Bulk RNA-Seq of Tlx3 WT vs KO samples
# Budhaditya Basu

library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(org.Mm.eg.db)
library(paletteer)
library(ComplexHeatmap)
# Load diff expressed gene Posterior cerebellum

Post <- read.csv("Post_WT_vs_KO.csv",
                 header = T, sep = "\t")%>%
  dplyr::select(Gene_name, logFC, adj.P.Val)

gene.list <- Post %>%
  arrange(adj.P.Val)%>%
  dplyr::mutate(Score = -log10(adj.P.Val + 0.0001)* sign(logFC))%>%
  dplyr::select(Gene_name, Score)

# Make the rank file
ranks <- deframe(gene.list)
head(ranks)

# Set decreasing order
geneList = sort(ranks, decreasing = TRUE)
head(geneList)
#===============================================================
# Perform GSEA
#===============================================================
m_ont <- msigdbr(species = "Mus musculus", 
                 category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)

head(m_ont)

em2 <- GSEA(geneList, TERM2GENE = m_ont)
head(em2)
gsea_result <- em2@result

# Save the GSEA result
write.csv(gsea_result, 
          file = "Post_GSEA_cell_type.csv", 
          row.names = F)

# Save the GSEA Plot
png(filename = "Granule_neuron_post_cerebellum.png", 
     width = 8, height = 6, units = "in",
     res = 600)
gseaplot2(em2, geneSetID = "DESCARTES_FETAL_CEREBELLUM_GRANULE_NEURONS", 
          title = em2$Description["DESCARTES_FETAL_CEREBELLUM_GRANULE_NEURONS"])+
  ggtitle("DESCARTES_FETAL_CEREBELLUM_GRANULE_NEURONS")
dev.off()

# DNA replication
png(filename = "DNA_replication_post_cerebellum.png", 
    width = 8, height = 6, units = "in",
    res = 600)
gseaplot2(em2, geneSetID = "GOBP_DNA_REPLICATION_CHECKPOINT_SIGNALING", 
          title = em2$Description["GOBP_DNA_REPLICATION_CHECKPOINT_SIGNALING"])+
  ggtitle("GOBP_DNA_REPLICATION_CHECKPOINT_SIGNALING")
dev.off()

png(filename = "WNT_signalling_post_cerebellum.png", 
    width = 8, height = 6, units = "in",
    res = 600)
gseaplot2(em2, geneSetID = "GOBP_NEGATIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY", 
          title = em2$Description["GOBP_NEGATIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY"])+
  ggtitle("GOBP_NEGATIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY")
dev.off()

# Matrix
png(filename = "Extracellular_matrix_post_cerebellum.png", 
    width = 8, height = 6, units = "in",
    res = 600)
gseaplot2(em2, geneSetID = "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", 
          title = em2$Description["GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX"])+
  ggtitle("GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX")
dev.off()

# Divergent Plot

Ont <- read.csv("Post_GSEA_Ontology.csv",
                header = T)
cell_type <- read.csv("Post_GSEA_cell_type.csv",
                      header = T)

Ont_top <- Ont %>%
  clusterProfiler::slice(3,7,8,24,26,67,35,36)%>% 
  mutate(gene_set = "Gene Ontology")

cell_top <- cell_type %>%
  clusterProfiler::slice(121,5,7,9,11)%>%
  mutate(gene_set = "Cell Type Signature")

df <- rbind(Ont_top, cell_top)

p <- ggplot(df, aes(x = reorder(ID, NES), y=NES))+
  geom_bar(stat='identity',  aes(fill=p.adjust), width=.8, color ="black",
           size = .1)+
  coord_flip()+
  theme_bw()+
  facet_grid(gene_set ~.,space="free", scales="free")+
  labs(x = "", y = "Normalized Enrichment Score", fill = "Adjusted P value")+
  paletteer::scale_fill_paletteer_c("grDevices::Blue-Yellow 3",
                                    direction = 1)+
  theme(axis.text.x =element_text(size=8, colour = "black"),
        axis.text.y = element_text(size=7),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        legend.title=element_text(size=6),
        legend.text=element_text(size=5),
        legend.key.width=unit(0.25,"cm"),
        legend.key.height=unit(0.25,"cm"),
        legend.position = c(0.7,0.1))+
  theme(strip.text.y = element_text(size = 8))
p

ggsave(filename = "Tlx3_KO_divergent.png",
       p,
       height = 5,
       width = 6,
       units = "in",
       dpi = 600)

# Heatmap
#========================================================================
# Load the data
data <- read.csv("Cerebellum_genes_final_normalized_scores.csv",
                 header = T)

head(data)

data <- data %>%
  separate(ID, into = c("ENS_ID", "Gene"), 
           sep = "_", remove = TRUE)

# Load diff expressed gene Anterior cerebellum

Post <- read.csv("Post_WT_vs_KO.csv",
                 header = T, sep = "\t")%>%
  select(Gene_name, logFC, adj.P.Val)

de_genes <- Post %>%
  filter(adj.P.Val < 0.05)

# Filter DE genes for heatmap

df <- data %>%
  filter(data$Gene %in% de_genes$Gene_name)%>%
  select(Gene, Ant_Cere_WT_N1, Ant_Cere_WT_N2,
         Ant_Cere_KO_N1, Ant_Cere_KO_N2,
         Post_Cere_WT_N1, Post_Cere_WT_N2,
         Post_Cere_KO_N1, Post_Cere_KO_N2)%>%
  column_to_rownames(var = "Gene")



mat <- as.matrix(df[,c(5:8)])
# Scale the expression value
scaled_mat <- t(scale(t(mat)))

my_palette <- paletteer::paletteer_c("grDevices::Blue-Yellow 3", 30)

ha = HeatmapAnnotation(df = data.frame(Samples = c("Control","Control",
                                                   "cKO","cKO")),
                       col = list(Samples = c("Control" = "#6FCEB4",
                                              "cKO" = "#DC0085")),
                       annotation_legend_param = list(title = "Post. Cerebellum",
                                                      title_gp = gpar(fontsize = 10,
                                                                      fontface = "bold")),
                       annotation_name_side = "left",
                       simple_anno_size = unit(2, "mm"))

ht1 <- Heatmap(scaled_mat,cluster_columns = T,
               width = unit(3, "cm"), 
               height = unit(8, "cm"),
               show_column_names = F,
               show_row_names = F,
               col = my_palette,
               bottom_annotation = ha,
               name = "Row Z score",
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10,
                                 fontface = "bold")
               ))

ht1

ann <- df %>%
  filter(rownames(df) %in% c("Grik2",
                             "Bmp5", "Ntf3", "Egr1",
                             "Ptpru", "Col6a2",
                             "Fbln2", "Col8a1",
                             "Cdc6", "Orc1"))

vrn <- rownames(df) %in% rownames(ann)

ht2 <- ht1 + 
  rowAnnotation(link = anno_mark(at = which(vrn), 
                                 labels = row.names(scaled_mat)[vrn],
                                 labels_gp = gpar(fontsize = 10), 
                                 padding = unit(1, "mm")))
ht2

png("Posterior/Post.Cerebellum_heatmap.png",
    height = 5, width = 6, units = "in",
    res = 600)
ht2
dev.off()