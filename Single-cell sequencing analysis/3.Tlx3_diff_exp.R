# Single cell RNA-Seq of Tlx3 WT vs KO samples
# Budhaditya Basu

library(Seurat)
library(tidyverse)
library(scCustomize)
library(SCP)
library(ggthemes)
library(RColorBrewer)
library(rstatix)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggrepel)
library(fgsea)
library(presto)
#=========================================================
tlx3.combined <- readRDS("sc_PN5_tlx3.combined_annotated.rds")
tlx3.combined@meta.data
DefaultAssay(tlx3.combined)
# Calculate the proportion of each cell type comparing the conditions
# How many cells are in each cluster
table(Idents(tlx3.combined))

# How many cells are in each replicate?
table(tlx3.combined$orig.ident)

# What proportion of cells are in each cluster?
prop.table(table(Idents(tlx3.combined)))

# How does cluster membership vary by condition?
table(Idents(tlx3.combined), tlx3.combined$orig.ident)


# Proportion of each cluster by condition
prop.table(table(Idents(tlx3.combined), tlx3.combined$orig.ident), margin = 2)

# Fisher's test to evaluate the power of cell composition analysis

test <- rstatix::row_wise_fisher_test(as.matrix(table(Idents(tlx3.combined), 
                                             tlx3.combined$orig.ident)),
                             p.adjust.method = "BH")

test
write.csv(test, file = "Fisher_test_for_cell_composition_analysis.csv",
          row.names = F)
# Keep the proportion in a data frame for plotting
df <- as.data.frame(prop.table(table(Idents(tlx3.combined), 
                                     tlx3.combined$orig.ident), margin = 2))

head(df)

plot <- ggplot(df, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_tufte() +
  geom_bar(position = "fill", width = 0.5, stat = "identity") +
  xlab("") +
  ylab("Proportion")+
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90))+
  scale_y_continuous(labels = scales::percent)+
  geom_text(aes(label = paste0(round(Freq*100, 2),"%"), y = Freq, x = Var2), 
            position = position_fill(vjust = 0.5), size = 2)

plot
ggsave(filename = "Plots/CelltypeProportion.jpg",
       plot,
       height = 6,
       width = 6,
       units = "in",
       dpi = 600)

# Granule Neuron Progenitor
#==========================================================================
#Differential Gene Expression Analysis of a cluster between two condition
markers <- FindMarkers(object = tlx3.combined, ident.1 = "KO", 
                       ident.2 = "Control", 
                       subset.ident = "Granule Neuron Prog", 
                       group.by = "orig.ident",
                       verbose = FALSE)
head(markers)
#Convert rownames to column names
df <- rownames_to_column(markers, var = "gene")
df
#Save the data frame as CSV file
write.csv(df, file = "Diff_exp_Granule_Neuron_Prog_WT_vs_KO.csv", row.names = F)

#Data wrangling
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
p <- vol+ 
  geom_point(size = 1, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#4a235a", 
                                Decreased = "#9a7d0a", 
                                nonsignificant = "gray80")) +
  theme_bw() + # change overall theme
  theme(legend.position = "none") + # change the legend
  xlab(expression(log[2]("KO" / "Control"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  scale_y_continuous(trans = "log1p")+ # Scale yaxis due to large p-values
  geom_hline(yintercept = 1.3, colour = "red", linetype="dashed") # Add p-adj value cutoff line
  #geom_vline(aes(xintercept=0.58), colour="gray60", linetype="dashed")+
  #geom_vline(aes(xintercept=-0.58), colour="gray60", linetype="dashed")
  #xlim(-2.5, 2.5)

# Base vocano plot
p



p2 <- p + 
  geom_text_repel(data=head(data,10), aes(label=gene),
                  box.padding = unit(0.5, "lines"),hjust= 0.30,
                  segment.color = 'black', max.overlaps = Inf,
                  colour = 'black', size = 3)+
  geom_text_repel(data = data %>%
                    filter(gene %in% c("Chrna3", 
                                       "Car8",
                                       "Sp5")), 
                  aes(label = gene, x = lfc, y = pval),
                  box.padding = unit(.5, "lines"),hjust= 0.30,
                  size = 3,
                  segment.color = 'black',
                  colour = 'black')+
  ggtitle("DGE GNP")
p2





ggsave(
  "Plots/volcano_GNP_KO_vs_WT.jpg",
  p2,
  width = 3.50,
  height = 4.00,
  dpi = 600
)



#Filter the genes which are upregulated in KO
geneList <- data %>% dplyr::filter(lfc >0 & pval > 1.3) 

ego <- enrichGO(gene          = geneList$gene,
                OrgDb         = org.Mm.eg.db, 
                ont           = "ALL", 
                pAdjustMethod = "fdr", 
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                keyType = "SYMBOL",
                readable      = TRUE)
head(ego)
head(ego@result)
write.csv(ego@result, file = "GO_GCP_clusterprofiler.csv", row.names = F)
#Plot 
tiff(filename = "GO_Granule_neuron_prog.tiff", width = 6, height = 6, units = "in",
     res = 600)
dotplot(ego)
dev.off()


#Filter the genes which are downregulated in KO
geneList_down <- data %>% dplyr::filter(lfc <0 & pval > 1.3) 

ego_down <- enrichGO(gene          = geneList_down$gene,
                OrgDb         = org.Mm.eg.db, 
                ont           = "ALL", 
                pAdjustMethod = "fdr", 
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                keyType = "SYMBOL", 
                readable      = TRUE)
head(ego)
head(ego@result)
write.csv(ego_down@result, file = "GO_GCP_downregulated.csv", row.names = F)

#Plot 
tiff(filename = "GO_downregulated_Granule_neuron_prog.tiff", width = 6, height = 6, units = "in",
     res = 600)
dotplot(ego_down)
dev.off()

#------------------------------------------------------------------------------

#Overrepresentation Analysis (ORA) using clusterprofiler by enricher()
msigdbr_species()
#Download MSigDb ontology gene sets
m_ont <- msigdbr(species = "Mus musculus", category = "C5") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_ont)

#MSigDb over-presentaton analysis
em <- enricher(gene = geneList$gene, TERM2GENE=m_ont)
head(em)

barplot(em)
dotplot(em)
upsetplot(em)

#Create upsetplot for overrepresentation analysis
tiff(filename = "upset.tiff", width = 16.00, height = 6.00,
     units = "in", res = 600)
upsetplot(em, 8)
dev.off()

# Granule Neuron
#=========================================================================
#Differential Gene Expression Analysis of a cluster between two condition
markers <- FindMarkers(object = tlx3.combined, ident.1 = "KO", 
                       ident.2 = "Control", 
                       subset.ident = "Granule Neuron", 
                       group.by = "orig.ident",
                       test.use = "wilcox",
                       verbose = FALSE)
head(markers)
#Convert rownames to column names
df <- rownames_to_column(markers, var = "gene")
df
#Save the data frame as CSV file
write.csv(df, file = "Diff_exp_Granule_Neuron_WT_vs_KO.csv", row.names = F)

#Data wrangling
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
p <- vol+ 
  geom_point(size = 1, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#4a235a", 
                                Decreased = "#9a7d0a", 
                                nonsignificant = "gray80")) +
  theme_bw() + # change overall theme
  theme(legend.position = "none") + # change the legend
  xlab(expression(log[2]("KO" / "Control"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  scale_y_continuous(trans = "log1p")+ # Scale yaxis due to large p-values
  geom_hline(yintercept = 1.3, colour = "red", linetype="dashed")  # Add p-adj value cutoff line
#geom_vline(aes(xintercept=0.58), colour="gray60", linetype="dashed")+
#geom_vline(aes(xintercept=-0.58), colour="gray60", linetype="dashed")+
#xlim(-2.5, 2.5)

#Base vocano plot
p

p2 <- p+ 
  geom_text_repel(data=head(data,10), aes(label=gene),
                  box.padding = unit(.5, "lines"),hjust= 0.30,
                  segment.color = 'black', max.overlaps = Inf,
                  colour = 'black', size = 3)+
  geom_text_repel(data = data %>%
                    filter(gene %in% c("Chrna3","Cre",
                                       "Nkd1", "Tle1",
                                       "Egr1", "Axin2", "Mllt3")), 
                  aes(label = gene, x = lfc, y = pval),
                  box.padding = unit(.5, "lines"),hjust= 0.30,
                  segment.color = 'black',
                  colour = 'black', size = 3)+
  ggtitle("DGE GN")
p2


ggsave(
  "Plots/volcano_Granule_Neuron_KO_vs_WT.tiff",
  p2,
  width = 3.50,
  height = 4.00,
  dpi = 600
)

# Filter the genes which are upregulated in KO
geneList <- data %>% dplyr::filter(lfc >0 & pval > 1.3) 

ego <- enrichGO(gene          = geneList$gene,
                OrgDb         = org.Mm.eg.db, 
                ont           = "ALL", 
                pAdjustMethod = "BH", 
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                keyType = "SYMBOL", 
                readable      = TRUE)
head(ego)
head(ego@result)
ego.result <- ego@result
write.csv(ego@result, file = "GO_Granule_neuron_clusterprofiler.csv", row.names = F)

#Plot 
tiff(filename = "GO_Granule_neuron.tiff", width = 10, height = 8, units = "in",
     res = 600)
dotplot(ego, showCategory = 20, font.size = 10)
dev.off()

# MSigDb over-presentaton analysis
em <- enricher(gene = geneList$gene, TERM2GENE = m_ont, pAdjustMethod = "fdr")
head(em)


#Purkinjee Neuron
#=========================================================================

#Differential Gene Expression Analysis of a cluster between two condition
markers <- FindMarkers(object = tlx3.combined, ident.1 = "KO", 
                       ident.2 = "Control", 
                       subset.ident = "Purkinje Cells", 
                       group.by = "orig.ident",
                       test.use = "wilcox",
                       verbose = FALSE)
head(markers)
#Convert rownames to column names
df <- rownames_to_column(markers, var = "gene")
df
#Save the data frame as CSV file
write.csv(df, file = "Diff_exp_Purkinji_WT_vs_KO.csv", row.names = F)

#Data wrangling
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
p <- vol+ 
  geom_point(size = 3, alpha = 0.5, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#4a235a", 
                                Decreased = "#9a7d0a", 
                                nonsignificant = "gray80")) +
  theme_classic() + # change overall theme
  theme(legend.position = "none") + # change the legend
  xlab(expression(log[2]("KO" / "Control"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  scale_y_continuous(trans = "log1p")+ # Scale yaxis due to large p-values
  geom_hline(yintercept = 1.3, colour = "red", linetype="dashed")  # Add p-adj value cutoff line
#geom_vline(aes(xintercept=0.58), colour="gray60", linetype="dashed")+
#geom_vline(aes(xintercept=-0.58), colour="gray60", linetype="dashed")+
#xlim(-2.5, 2.5)

#Base vocano plot
p

p2 <- p+ geom_text_repel(data=head(data,10), aes(label=gene),
                         box.padding = unit(.2, "lines"),hjust= 0.30,
                         segment.color = 'black', max.overlaps = Inf,
                         colour = 'black')+
  geom_text_repel(data = data %>%
                    filter(gene %in% c("Chrna3")), 
                  aes(label = gene, x = lfc, y = pval),
                  box.padding = unit(.2, "lines"),hjust= 0.30,
                  segment.color = 'black',
                  colour = 'black')
p2


ggsave(
  "volcano_Granule_Neuron_KO_vs_WT.tiff",
  p2,
  width = 6.00,
  height = 4.00,
  dpi = 600
)


# Filter the genes which are upregulated in KO
geneList <- data %>% dplyr::filter(lfc >0) 

ego <- enrichGO(gene          = geneList$gene,
                OrgDb         = org.Mm.eg.db, 
                ont           = "ALL", 
                pAdjustMethod = "BH", 
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                keyType = "SYMBOL", 
                readable      = TRUE)
head(ego)

# Plot 
tiff(filename = "GO_Granule_neuron.tiff", width = 6, height = 6, units = "in",
     res = 600)
dotplot(ego, showCategory = 10)
dev.off()

# MSigDb over-presentaton analysis
em <- enricher(gene = geneList$gene, TERM2GENE = m_ont, pAdjustMethod = "fdr")
head(em)


#Perform GSEA analysis on differential Gene Expression Analysis
#==============================================================================

markers <- FindMarkers(object = tlx3.combined, ident.1 = "KO", 
                       ident.2 = "Control", 
                       subset.ident = "Granule Neuron", 
                       group.by = "orig.ident",
                       test.use = "wilcox",
                       verbose = FALSE)
head(markers)
# Convert rownames to column names
df <- rownames_to_column(markers, var = "gene")

gene.list <- df %>%
  arrange(p_val_adj)%>%
  mutate(score = -log10(p_val_adj)*sign(avg_log2FC))%>%
  dplyr::select(gene, score)

#Make the rank file
ranks <- deframe(gene.list)

head(ranks)

#we need to have the annotated gene set first. 

msigdbr_show_species()

m_df <- msigdbr(species = "Mus musculus", category = "C5")

#filtering and formating gene sets
fgsea_sets <- m_df %>% split(x = .$gene_symbol, 
                             f = .$gs_name)
#Perform GSEA
fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

#data wrangling
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

head(fgseaResTidy)


# GSEA plot
tiff(filename = "GSEA_NSC_vs_RGC_GO_Neuron_diff.tiff", height = 3, width = 6,
     units = "in", res = 200)
plotEnrichment(fgsea_sets[["GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION"]],
               ranks) +
  labs(title="Cell Morphogenesis in NEURON_DIFFERENTIATION")
dev.off()

gcp_upregulated <- read.csv("Diff_exp_Granule_Neuron_Prog_WT_vs_KO.csv")%>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0)%>%
  select(gene)
write.csv(gcp_upregulated, file = "gcp_upregulated_genes.csv")
head(gcp_upregulated)

gcp_downregulated <- read.csv("Diff_exp_Granule_Neuron_Prog_WT_vs_KO.csv")%>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) 

gn_upregulated <- read.csv("Diff_exp_Granule_Neuron_WT_vs_KO.csv")%>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0)%>%
  select(gene)

write.csv(gn_upregulated, file = "gn_upregulated_genes.csv")

read.csv("Diff_exp_Granule_Neuron_WT_vs_KO.csv")%>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
  nrow()

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
