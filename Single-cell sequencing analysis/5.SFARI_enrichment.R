# Single cell RNA-Seq of Tlx3 WT vs KO samples
# Budhaditya Basu

install.packages("BiocManager")
BiocManager::install("AnnotationDbi")

#To install the geneLists package in R, use the devtools package:
  
# Install from github
devtools::install_github("soderling-lab/geneLists")
devtools::install_github("montilab/hypeR")
#The gene mapping function getIDs uses organism specific mapping data. 
  
BiocManager::install("org.Mm.eg.db")

#==========================================================================
## using geneLists to perform GSEA
library(geneLists)
library(tidyverse)
library(org.Mm.eg.db)
library(biomaRt)
library(ggrepel)
library(ggsci)
# to see all available gene lists:
# geneLists()
# 
# msigdb_info()
# KEGG <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory = "CP:KEGG")
# GO.BP <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory = "BP")
# GO.CC <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory = "CC")
# GO.MF <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory = "MF")

# load the SFARI ASD-associated genes
data("sfariGene", package = "geneLists")
data("sfariAnimal", package = "geneLists")
# a geneList is just a named list of Entrez IDs
str(sfariGene)

# combine SFARI gene lists
sfari <- unique(c(sfariGene[["ASD"]], sfariAnimal[["ASD"]]))

#Gene set enrichment and over-representation
#===============================================================================
#Identify DEGs between Control and KO
#=========================================================
#Import DEGs of three clusters (GN, GNP and Purkinji cells)
gn <- read.csv(file = "Diff_exp_Granule_Neuron_WT_vs_KO.csv", header = T)%>%
  filter(p_val_adj < 0.05)
gnp <- read.csv(file = "Diff_exp_Granule_Neuron_Prog_WT_vs_KO.csv", header = T)%>%
  filter(p_val_adj < 0.05)
pn <- read.csv(file = "Diff_exp_Purkinji_WT_vs_KO.csv", header = T)%>%
  filter(p_val_adj < 0.05)

#Combine the data
data <- rbind(gn, gnp, pn)%>%
  distinct()

#Convert the gene symbol to Entrez id
entrez_id <- mapIds(org.Mm.eg.db, data$gene, 'ENTREZID', 'SYMBOL')

#===============================================================================
#Import all mouse entez_id as background for hypergeometric test
# Connect to the Ensembl database using the "mmusculus_gene_ensembl" dataset
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get all mouse gene Entrez IDs using the "entrezgene" attribute
mouse_entrez_ids <- getBM(attributes = c("entrezgene_id"),
                          mart = ensembl)

#===============================================================================
# Perform hypergeometric test for enrichment of SFARI genes in Tlx3 KO data
hyp_obj <- geneLists::hyperTest(sfari, entrez_id, 
                                background = mouse_entrez_ids$entrezgene_id)
head(hyp_obj)
write.csv(hyp_obj, file = "Hypergeometric_test_SFARI.csv")

# convert entrez IDs to mouse gene symbol using getIDs
sfari_genes <- geneLists::getIDs(sfari, "entrez", "symbol", "mouse")

overlap_genes <- sfari_genes[which(sfari_genes %in% data$gene)]%>%
  stringr::str_to_upper()


#Import SFARI genes
sfari_list <- read.csv("SFARI-Gene.csv", header = T)

#Collect overlapping genes data
collect_data <- sfari_list %>%
  filter(gene.symbol %in% overlap_genes)

write.csv(collect_data, file = "Tlx3_data_OverlapSFARI_database.csv",
          row.names = F)

collect_data$gene.symbol <- str_to_title(collect_data$gene.symbol)
#Make a plot

p <- ggplot(collect_data, aes(x = status, y = number.of.reports,
                                  fill = factor(syndromic)))+
  geom_point(shape = 21,
             aes(size = factor(gene.score)))+
  theme_bw()+
  labs(x = "Genes",
       y = "Number of Reports")+
  scale_x_continuous(breaks = c(8, 12))+
  geom_text_repel(aes(label = gene.symbol),
                  size = 3,
                  box.padding = unit(.5, "lines"),hjust= 0.30)+
  scale_fill_jco(name = "Syndromic", 
                 labels = c("Yes", 
                            "No"))+
  scale_size_discrete(name = "Gene Score",
                      labels = c("1 (High Confidence)",
                                 "2 (Strong Candidate)",
                                 "Missing Value"))+
  geom_label(
    label= "Hypergeometric test\nFold enrichment = 2.08\nP-Value=0.023", 
    x=9.03,
    y=27,
    label.size = 0.01,
    show.legend = F)+
  ggtitle("Overlap of DEGs with SFARI ASD database")

p
ggsave(filename = "Plots/SFARI_ASD_overlap.jpg",
       p,
       height = 7,
       width = 7,
       units = 'in',
       dpi = 600)
