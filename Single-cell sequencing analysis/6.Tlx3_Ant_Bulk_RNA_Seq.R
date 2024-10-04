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

# Load diff expressed gene Anterior cerebellum

Ant <- read.csv("Ant_WT_vs_KO.csv",
                 header = T, sep = "\t")%>%
  dplyr::select(Gene_name, logFC, adj.P.Val)

gene.list <- Ant %>%
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
                 category = "C8") %>% 
  dplyr::select(gs_name, gene_symbol)

head(m_ont)

em2 <- GSEA(geneList, TERM2GENE = m_ont)
head(em2)
gsea_result <- em2@result

# Save the GSEA result
write.csv(gsea_result, 
          file = "Ant_GSEA_cell_type.csv", 
          row.names = F)

# Save the GSEA Plot
png(filename = "GOBP_CYTOPLASMIC_TRANSLATION.png", 
    width = 8, height = 6, units = "in",
    res = 600)
gseaplot2(em2, geneSetID = "GOBP_CYTOPLASMIC_TRANSLATION", 
          title = em2$Description["GOBP_CYTOPLASMIC_TRANSLATION"])+
  ggtitle("GOBP_CYTOPLASMIC_TRANSLATION")
dev.off()

png(filename = "GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY.png", 
    width = 8, height = 6, units = "in",
    res = 600)
gseaplot2(em2, geneSetID = "GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY", 
          title = em2$Description["GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY"])+
  ggtitle("GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY")
dev.off()

# Divergent Plot

Ont <- read.csv("Ant_GSEA_Ontology.csv",
                header = T)
cell_type <- read.csv("Ant_GSEA_cell_type.csv",
                      header = T)

Ont_top <- Ont %>%
  clusterProfiler::slice(1,5,7,10,13,15,17,25,160)%>% 
  mutate(gene_set = "Gene Ontology")

cell_top <- cell_type %>%
  clusterProfiler::slice(1,3,6,7, 70, 108, 124)%>%
  mutate(gene_set = "Cell Type Signature")

df <- rbind(Ont_top, cell_top)

p <- ggplot(df, aes(x = reorder(ID, NES), y=NES))+
  geom_bar(stat='identity',  aes(fill=p.adjust), width=.8, color ="black",
           size = .1)+
  coord_flip()+
  theme_bw()+
  facet_grid(gene_set ~.,space="free", scales="free")+
  labs(x = "", y = "Normalized Enrichment Score", fill = "Adjusted P value")+
  paletteer::scale_fill_paletteer_c("grDevices::Cyan-Magenta",
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

ggsave(filename = "Tlx3_KO_Ant_divergent.png",
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

Ant <- read.csv("Ant_WT_vs_KO.csv",
                header = T, sep = "\t")%>%
  select(Gene_name, logFC, adj.P.Val)

de_genes <- Ant %>%
  filter(adj.P.Val < 0.05)

# Filter DE genes for heatmap

df <- data %>%
  filter(data$Gene %in% de_genes$Gene_name)%>%
  select(Gene, Ant_Cere_WT_N1, Ant_Cere_WT_N2,
         Ant_Cere_KO_N1, Ant_Cere_KO_N2,
         Post_Cere_WT_N1, Post_Cere_WT_N2,
         Post_Cere_KO_N1, Post_Cere_KO_N2)%>%
  column_to_rownames(var = "Gene")



mat <- as.matrix(df[,c(1:4)])
# Scale the expression value
scaled_mat <- t(scale(t(mat)))

my_palette <- paletteer::paletteer_c("grDevices::Cyan-Magenta", 30)

ha = HeatmapAnnotation(df = data.frame(Samples = c("Control","Control",
                                                   "cKO","cKO")),
                       col = list(Samples = c("Control" = "khaki3",
                                              "cKO" = "darkmagenta")),
                       annotation_legend_param = list(title = "Ant. Cerebellum",
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
  filter(rownames(df) %in% c("Rpl38", "Rpl36a",
                             "Kcna1", "Lrrc26",
                             "Cacna2d1", "Ano3",
                             "Hmgb1", "Hmgb2", "Stmn1",
                             "Cenpe", "Birc5", "Pcsk9"))
vrn <- rownames(df) %in% rownames(ann)

ht2 <- ht1 + 
  rowAnnotation(link = anno_mark(at = which(vrn), 
                                 labels = row.names(scaled_mat)[vrn],
                                 labels_gp = gpar(fontsize = 10), 
                                 padding = unit(1, "mm")))
ht2

png("Anterior/Ant.Cerebellum_heatmap.png",
    height = 5, width = 6, units = "in",
    res = 600)
ht2
dev.off()
