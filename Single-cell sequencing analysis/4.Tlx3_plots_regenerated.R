# Single cell RNA-Seq of Tlx3 WT vs KO samples
# Budhaditya Basu
library(Seurat)
library(tidyverse)
library(scCustomize)
library(SCP)
library(ggrepel)
library(ggthemes)
# =========================================================
tlx3.combined <- readRDS("sc_PN5_tlx3.combined_annotated.rds")
head(tlx3.combined@meta.data)
# Check Cre expression and Atoh1 expression
DefaultAssay(tlx3.combined)
# DefaultAssay(tlx3.combined) <- "RNA"

tlx3.combined$CellType <- Idents(tlx3.combined)

# Save split_dimplot using SCP
jpeg(filename = "Plots/Split_dimplot_SCP_Tlx3.jpg",
     height = 10, width = 15, units = 'in', res = 600)
CellDimPlot(
  srt = tlx3.combined, split.by = 'orig.ident', 
  group.by = 'CellType',
  reduction = "UMAP", theme_use = "theme_blank",
  pt.alpha = 2, pt.size = 2,
  legend.position = 'bottom')
dev.off()

# Save grouped_dimplot using SCP
jpeg(filename = "Plots/grouped_dimplot_SCP_Tlx3.jpg",
     height = 10, width = 12, units = 'in', res = 600)
CellDimPlot(
  srt = tlx3.combined,
  group.by = 'orig.ident',
  reduction = "UMAP", theme_use = "theme_blank",
  pt.alpha = 2, pt.size = 2,
  legend.position = 'right')
dev.off()

# Cell Cycle Scoring
# Load cell cycle scoring gene set
load("~/Documents/Single_cell_data/Mouse_cell_cycle_genes.RData")

# Score cells for cell cycle
tlx3.combined <- CellCycleScoring(tlx3.combined,
                                  g2m.features = g2m_genes,
                                  s.features = s_genes)

head(tlx3.combined@meta.data)

# Cell_cycle scoring
?CellDimPlot

jpeg(filename = "Plots/Cell_cycle_SCP_Tlx3.jpg",
     height = 10, width = 12, units = 'in', res = 600)
CellDimPlot(
srt = tlx3.combined, group.by = "CellType", stat.by = "Phase",
reduction = "UMAP", theme_use = "theme_blank", stat_plot_size = 0.05)
dev.off()


#========================================================
# GroupHeatmap
genes <- list("GNP" = c("Mfap4", "Atoh1", "Tlx3"),
              "GN" = c("Neurod1", "Cntn2"),
              "GABAergic" = c("Pax2"),
              "Purkinji" = c("Calb1", "Car8"),
              "Bergmann Glia" = c("Hopx", "Gdf10"),
              "Oligodendrocyte" = c("Sox10", "Olig1"),
              "Endo Precursor" = c("Lgals1", "Col3a1"),
              "Erythrocytes" = c("Hbb-bt", "Hba-a2"),
              "Endothelial" = c("Cldn5", "Kdr"))
ht <- GroupHeatmap(
  srt = tlx3.combined,
  features = c(
    "Mfap4", "Atoh1", "Tlx3", # GNP
    "Neurod1", "Cntn2", # GN
    "Pax2", # GABAergic
    "Calb1", "Car8", # Purkinji
    "Hopx", "Gdf10", # Burgmann
    "Sox10", "Olig1", # Oligodendrocyte
    "Lgals1", "Col3a1", #Endo precursor
    "Hbb-bt", "Hba-a2", #Erythrocyte
    "Cldn5", "Kdr", #Endothelial
    "Chrna3"
  ),
  group.by = c("CellType"),
  split.by = "orig.ident",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score", "Cre"),
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
print(ht$plot)

# Save this plot
jpeg(filename = "Plots/Markers_SCP_Tlx3.jpg",
     height = 8, width = 15, units = 'in', res = 600)
ht$plot
dev.off()

#============================================================================
#Differential expression analysis

# Check all disturbed marker in KO samples
markers <- FindMarkers(object = tlx3.combined, ident.1 = "KO", 
                       ident.2 = "Control", 
                       group.by = "orig.ident",
                       verbose = FALSE)

# Convert rownames to column names
df <- rownames_to_column(markers, var = "gene")
# Save the data frame as CSV file
write.csv(df, file = "Diff_exp_allcluster_WT_vs_KO.csv", row.names = F)

#===================================================================
# DETest in SCP
tlx3.combined <- RunDEtest(srt = tlx3.combined,
                           group_by = "CellType")

df <- tlx3.combined@tools$DEtest_CellType$AllMarkers_wilcox
markers <- filter(df, p_val_adj < 0.05 & avg_log2FC > 1)

table(markers$group1)

# Annotate features with transcription factors and surface proteins
tlx3.combined <- AnnotateFeatures(tlx3.combined,
                                  species = "Mus_musculus",
                                  db = c("TF", "SP"))
ht <- FeatureHeatmap(
  srt = tlx3.combined, group.by = "CellType",
  features = markers$gene, feature_split = markers$group1,
  species = "Mus_musculus",
  db = c("GO_BP", "KEGG", "WikiPathway"),
  anno_terms = TRUE,
  feature_annotation = c("TF", "SP"),
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 10, width = 6)

print(ht$plot)

jpeg(filename = 'Plots/FeatureHeatmap_between_clusters_new.jpg',
     height = 12, width = 25, units = 'in', res = 600)
ht$plot
dev.off()

#=======================================================================
# Expression difference

jpeg(filename = 'Plots/Featurestatplot_split.jpg',
     height = 10, width = 10, units = 'in', res = 600)
FeatureStatPlot(
  srt = tlx3.combined, group.by = "CellType", bg.by = "CellType", 
  split.by = 'orig.ident',
  stat.by = c("Cre", "Chrna3", "Calb1", "Car8", "Sp5",
              "Ctnnb1", "Rora", "Atp1a3"), add_box = TRUE,
  comparisons = list(
    c("Control", "KO")),
  stack = TRUE, xlab = NULL)
dev.off()

?FeatureStatPlot

# ASD gene list
jpeg(filename = 'Plots/Featurestatplot_ASD_genes_split.jpg',
     height = 10, width = 10, units = 'in', res = 600)
FeatureStatPlot(
  srt = tlx3.combined, group.by = "CellType", bg.by = "CellType", 
  split.by = 'orig.ident',
  stat.by = c("Ctnnb1", "Rora", "Atp1a3",
              "Dip2a", "Cadm1", "Snap25",
              "Ptbp2",
              "Zswim6", "Dvl3", "Rheb"), add_box = TRUE,
  comparisons = list(
    c("Control", "KO")),
  stack = TRUE, xlab = NULL)
dev.off()


# Differential gene expression
#========================================================================
markers <- FindMarkers(object = tlx3.combined, ident.1 = "KO", 
                       ident.2 = "Control", 
                       subset.ident = "Purkinje Cells", 
                       group.by = "orig.ident",
                       test.use = "wilcox",
                       verbose = FALSE)
head(markers)
# Convert rownames to column names
df <- rownames_to_column(markers, var = "gene")
df
# Save the data frame as CSV file
write.csv(df, file = "Diff_exp_Purkinji_WT_vs_KO.csv", row.names = F)

# Data wrangling
data <- data.frame(gene = row.names(markers),
                   pval = -log10(markers$p_val_adj), 
                   lfc = markers$avg_log2FC)


data <- mutate(data, color = case_when(data$lfc > 0 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 0 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))

head(data)

# Make a basic ggplot2 object with x-y values
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))


# Add ggplot2 layers (COLORS : khaki =#9a7d0a, Purple = #4a235a )
p <- vol + 
  geom_point(size = 1, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#4a235a", 
                                Decreased = "#9a7d0a", 
                                nonsignificant = "gray80")) +
  theme_bw() + # change overall theme
  theme(legend.position = "none") + # change the legend
  xlab(expression(log[2]("KO" / "Control"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  #scale_y_continuous(trans = "log1p")+ # Scale yaxis due to large p-values
  geom_hline(yintercept = 1.3, colour = "red", linetype="dashed")+ # Add p-adj value cutoff line
  
  geom_vline(aes(xintercept=0.58), colour="gray60", linetype="dashed")+
  geom_vline(aes(xintercept=-0.58), colour="gray60", linetype="dashed")
  #xlim(-2.5, 2.5)

# Base vocano plot
p

p2 <- p+ geom_text_repel(data=head(data,9), aes(label=gene),
                         box.padding = unit(.2, "lines"),hjust= 0.30,
                         segment.color = 'black', max.overlaps = Inf,
                         colour = 'black', size = 3)+
  geom_text_repel(data = data %>%
                    filter(gene %in% c("Chrna3")), ##For mouse use title case
                  aes(label = gene, x = lfc, y = pval),
                  box.padding = unit(.2, "lines"),hjust= 0.30,
                  segment.color = 'black',
                  colour = 'black')+
  ggtitle("DGE Purkinji Cells")
p2


ggsave(
  "Plots/volcano_Purkinji_KO_vs_WT.jpg",
  p2,
  width = 3.50,
  height = 4.00,
  dpi = 600
)


# DAVID dotplot

gcp_david <- read.csv(file = "DAVID_gcp_upreg_new.txt",
                      sep = "\t", skip = 1)
gn_david <- read.csv(file = "DAVID_gn_upreg_new.txt", sep = "\t", skip = 1)

selected_gcp_david <- gcp_david[c(1,2,6,8,9), ]%>%
  dplyr::select(Term, Fold.Enrichment, PValue, Count, Genes)%>%
  mutate(Celltype = "GCP")

selected_gn_david <- gn_david[c(11,12,13,24,38,39, 55, 56), ]%>%
  dplyr::select(Term, Fold.Enrichment, PValue, Count, Genes)%>%
  mutate(Celltype = "GN")

df1 <- rbind(selected_gcp_david, selected_gn_david)
write.csv(df1, file = "DAVID_upregulated_dataset_GCP_GN.csv", row.names = F)

p <- ggplot(df1, aes(x = Celltype, y = Term, fill = as.numeric(PValue)))+
  geom_point(aes(size = as.numeric(Fold.Enrichment)),
             shape=21, colour = "black")+
  theme_bw()+  #Change theme to black and white
  scale_fill_gradient2(low="darkcyan",mid ="red", high="blue")+
  labs(size = "Fold Enrichment", x = "", y = "",
       fill = " P Value")
p

ggsave(filename = "DAVID_analysis.tiff",
       p,
       height = 5.00,
       width = 8.00,
       dpi = 600)


head(tlx3.combined@meta.data)
FeatureDimPlot(tlx3.combined, 
               features = c("Atoh1",
                            "Cre",
                            "Tlx3"),
               split.by = "orig.ident",
               reduction = "umap",
               theme_use = "theme_blank",
               nrow = 3)

png("Atoh1_dimplot.png", height = 4,
    width = 8, units = "in", res = 600)
FeaturePlot_scCustom(tlx3.combined, 
                     features = c("Atoh1"),
                     split.by = "orig.ident",
                     order = TRUE)& NoAxes()
dev.off()

png("Cre_dimplot.png", height = 4,
    width = 8, units = "in", res = 600)
FeaturePlot_scCustom(tlx3.combined, 
                     features = c("Cre"),
                     split.by = "orig.ident",
                     order = TRUE)& NoAxes()
dev.off()


png("Tlx3_dimplot.png", height = 4,
    width = 8, units = "in", res = 600)
FeaturePlot_scCustom(tlx3.combined, 
                     features = c("Tlx3"),
                     split.by = "orig.ident",
                     order = TRUE)& NoAxes()
dev.off()

#stacked violin

png("Stacked_violin_new.png",
    height = 6, width = 8, units = "in", res = 600)
Stacked_VlnPlot(tlx3.combined, 
                features = c("Atoh1",
                             "Cre",
                             "Tlx3"),
                split.by = "orig.ident",
                x_lab_rotate = TRUE,
                colors_use = c("khaki3", "darkmagenta"))
dev.off()


png("Stacked_violin.png",
    height = 6, width = 10, units = "in", res = 600)
FeatureStatPlot(
  srt = tlx3.combined, group.by = "CellType", 
  bg.by = "CellType", 
  split.by = 'orig.ident',
  stat.by = c("Atoh1", "Cre", "Tlx3"), 
  add_box = TRUE,
  stack = TRUE, xlab = "",
  palcolor = c("khaki3", "darkmagenta"))
dev.off()

#
png("Cadm1_dimplot.png", height = 4,
    width = 8, units = "in", res = 600)
FeaturePlot_scCustom(tlx3.combined, 
                     features = c("Cadm1"),
                     split.by = "orig.ident",
                     order = TRUE)& NoAxes()
dev.off()

png("Cdk2ap1_dimplot.png", height = 4,
    width = 8, units = "in", res = 600)
FeaturePlot_scCustom(tlx3.combined, 
                     features = c("Cdk2ap1"),
                     split.by = "orig.ident",
                     order = TRUE)& NoAxes()
dev.off()

png("Tgfb2_dimplot.png", height = 4,
    width = 8, units = "in", res = 600)
FeaturePlot_scCustom(tlx3.combined, 
                     features = c("Tgfb2"),
                     split.by = "orig.ident",
                     order = TRUE)& NoAxes()
dev.off()

png("Bbip1_dimplot.png", height = 4,
    width = 8, units = "in", res = 600)
FeaturePlot_scCustom(tlx3.combined, 
                     features = c("Bbip1"),
                     split.by = "orig.ident",
                     order = TRUE)& NoAxes()
dev.off()

png("Insc_dimplot.png", height = 4,
    width = 8, units = "in", res = 600)
FeaturePlot_scCustom(tlx3.combined, 
                     features = c("Insc"),
                     split.by = "orig.ident",
                     order = TRUE)& NoAxes()
dev.off()


png("Proliferation_new.png",
    height = 6, width = 8, units = "in", res = 600)
Stacked_VlnPlot(tlx3.combined, 
                features = c("Cadm1",
                             "Cdk2ap1",
                             "Tgfb2",
                             "Bbip1",
                             "Insc",
                             "Tlx3"),
                split.by = "orig.ident",
                x_lab_rotate = TRUE,
                colors_use = c("khaki3", 
                               "darkmagenta"))
dev.off()


