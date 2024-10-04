
#Find Tlx3 levels in different cell types of cerebellum at different dev stages
#Data is from : https://doi.org/10.1016/j.cub.2018.07.062
#Analysis by Rahul Jose
#
library(tidyverse)
library(gdata)
library(Seurat)

# Move and rename raw files -----------------------------------------------

#file lists
mtxlist <- list.files(path = "R:/Cerebellum_dev/raw/extract", pattern="\\.mtx.gz$")
alllist <- list.files(path = "R:/Cerebellum_dev/raw/extract") 
all_list2 <-  str_replace_all(alllist, pattern = "-", replacement = "_") 


listf <- data.frame("name"=all_list2, "name_old"=alllist)
 
#details from file name
listf <- listf |> 
  separate(col = name, sep = "_", into = c("id", "id2", "id3", "sample", "file"), remove = F)

listf <- mutate(listf, stage = str_sub(sample, end = -2))

# -------------------------------------------------------------------------

  
#create separate folders
for (x in 1:length(unique(listf$stage))) {
  
  dir.create(path = file.path("raw", "raw_sort", unique(listf$stage)[x]))
}

#copy respective files
  for(y in 1:length(listf$name)){
  file.copy(from = file.path("raw", "extract", listf$name_old[y]),
            to = file.path("raw", "raw_sort", listf$stage[y], listf$name[y]))
  }

#change names to just barcodes.tsv.gz, genes.tsv.gz and matrix.mtx.gz
for(y in 1:length(listf$name)){
  file.rename(from = file.path("raw", "raw_sort", listf$stage[y], listf$name[y]),
            to = file.path("raw", "raw_sort", listf$stage[y], listf$file[y]))
}

#change names from genes.tsv.gz to  features.tsv.gz
listf2 <- filter(listf, file=="genes.tsv.gz")
for(y in 1:length(listf2$name)){
  file.rename(from = file.path("raw", "raw_sort", listf2$stage[y], "genes.tsv.gz"),
              to = file.path("raw", "raw_sort", listf2$stage[y], "features.tsv.gz") )
}



# Read 10x files-------------------------------------------------------------------------


for(i in 1:length(listf2$stage)){
  temp1 <- Read10X(data.dir = file.path("raw", "raw_sort", listf2$stage[i]))
  temp2 <- CreateSeuratObject(temp1)
  mv(from = "temp2", to = tolower(listf2$stage[i]))
  rm(temp1)
}

#save(list = ls(),file = "Seurats.RData")

# Add cell labels as metadata ---------------------------------------------

seurats1 <- list("e10c"=e10c,"e10d"=e10d,"e11a"=e11a,"e11b"=e11b,"e12a"=e12a,"e12b"=e12b,
                 "e13a"=e13a,"e13b"=e13b,"e14a"=e14a,"e14b"=e14b,"e15a"=e15a,"e15b"=e15b,"e16c"=e16c,
                 "e16d"=e16d,"e17a"=e17a,"e17b"=e17b,"p0a"=p0a,"p0b"=p0b,"p10a"=p10a,"p10b"=p10b,
                 "p4a"=p4a,"p4b"=p4b,"p7a"=p7a,"p7b"=p7b)

names(seurats1[1])


listf3 <- filter(listf, file=="matrix.mtx.gz")

for (x in 1:length(listf3$stage)) {
  
  temp1 <- seurats1[[which(names(seurats1)== tolower(listf3$sample[x]))]]
  temp1 <- AddMetaData(temp1, col.name = "stage", metadata = listf3$stage[x])
  temp1 <- AddMetaData(temp1, col.name = "sample", metadata = listf3$sample[x])
  rm(list = c(tolower(listf3$sample[x])))
  mv(from = "temp1", to = tolower(listf3$sample[x]))
  #rm(temp1)
}
rm(seurats1)

cereb <- merge(x= e10c, y=c(e10d,e11a,e11b,e12a,e12b,
                            e13a,e13b,e14a,e14b,e15a,
                            e15b,e16c,e16d,e17a,e17b,
                            p0a,p0b,p10a,p10b,p4a,
                            p4b,p7a,p7b))
# saveRDS(object = cereb, file="cereb_merged.rds")
cereb1 <- readRDS(file="cereb_merged.rds")  


#QC -------------------------------------------------------------------------
#Reference : https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

genes <- rownames(cereb1@assays[["RNA"]]@features)
genes_df <- data.frame(gene=genes)

# Add number of genes per UMI for each cell to metadata
cereb1$log10GenesPerUMI <- log10(cereb1$nFeature_RNA) / log10(cereb1$nCount_RNA)

# Compute percent mito ratio
cereb1$mitoRatio <- PercentageFeatureSet(object = cereb1, pattern = "^mt-")
cereb1$mitoRatio <- cereb1@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- cereb1@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata |>
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
cereb1@meta.data <- metadata

# Create .RDS object to load at any time
saveRDS(cereb1, file ="seurat_cereb1_before-qc-metadata.rds")

########QC for
# Cell counts
# UMI counts per cell
# Genes detected per cell
# UMIs vs. genes detected
# Mitochondrial counts ratio
# Novelty
# Doublets

# Visualize the number of cell counts per sample
metadata |>
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("No. of Cells")
ggsave((filename = "before_qc_no_of_cells.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)


# Visualize the number UMIs/transcripts per cell

metadata |>
  ggplot(aes(color=sample, x=nUMI
             # , fill= sample
             )) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave((filename = "before_qc2_UMI_counts_per_cell.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

# Genes detected per cell
#Visualize the distribution of genes detected per cell via histogram
metadata |>
  ggplot(aes(color=sample, x=nGene
             # , fill= orig.ident
             )) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)
ggsave((filename = "before_qc3_genes_per_cell_histo.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)
# Visualize the distribution of genes detected per cell via boxplot
metadata |>
  ggplot(aes(x=sample, y=log10(nGene)
             # , fill=orig.ident
             )) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave((filename = "before_qc4_genes_per_cell_box.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata |>
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~stage)
ggsave((filename = "before_qc5_genes_umi_mito.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata  |>
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.2)
ggsave((filename = "before_qc6_mito_per_cell.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

#Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata |>
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave((filename = "before_qc7_genes_per_umi.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

# Filtering and QC-------------------------------------------------------------------------

##### Cell-level filtering
# nUMI > 500
# nGene > 250
# log10GenesPerUMI > 0.8
# mitoRatio < 0.2

cereb2 <- subset(x = cereb1, 
                 subset= (nUMI >= 500) & 
                   (nGene >= 250) & 
                   (log10GenesPerUMI > 0.75) &
                   (mitoRatio < 0.20))
metadata2 <- cereb2@meta.data

#POST-FILTERING QC
# Visualize the number of cell counts per sample
metadata2 |>
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("No. of Cells")
ggsave((filename = "after_qc_no_of_cells.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)


# Visualize the number UMIs/transcripts per cell

metadata2 |>
  ggplot(aes(color=sample, x=nUMI
             # , fill= sample
  )) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
ggsave((filename = "after_qc2_UMI_counts_per_cell.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

#Visualize the distribution of genes detected per cell via histogram
metadata2 |>
  ggplot(aes(color=sample, x=nGene
             # , fill= orig.ident
  )) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)
ggsave((filename = "after_qc3_genes_per_cell_histo.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)
# Visualize the distribution of genes detected per cell via boxplot
metadata2 |>
  ggplot(aes(x=sample, y=log10(nGene)
             # , fill=orig.ident
  )) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave((filename = "after_qc4_genes_per_cell_box.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata2 |>
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~stage)
ggsave((filename = "after_qc5_genes_umi_mito.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata2  |>
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.2)
ggsave((filename = "after_qc6_mito_per_cell.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)

#Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata2 |>
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
ggsave((filename = "after_qc7_genes_per_umi.jpeg"),
       plot=last_plot(), path = "plots/QC/", dpi = 320
)


# -------------------------------------------------------------------------
#Gene-level filtering
# Output a logical vector for every gene on whether the more than zero counts per cell

#join layers
cereb2_1 <- JoinLayers(object = cereb2)
# Extract counts
counts <- LayerData(object = cereb2_1, layer = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat2 <- CreateSeuratObject(filtered_counts, meta.data = cereb2_1@meta.data)
saveRDS(filtered_seurat2, file = "filtered_seurat2.rds")

#sample e12a is not good; excluding e12a
cerebx <- readRDS("filtered_seurat2.rds")
qc_passed <- toupper(c("e10c", "e10d","e11a","e11b",
#e12a,
"e12b","e13a","e13b","e14a","e14b","e15a","e15b","e16c","e16d","e17a","e17b","p0a","p0b","p10a","p10b","p4a","p4b","p7a","p7b"))
cerebxs <- subset(x = cerebx, subset = sample %in% qc_passed)
saveRDS(cerebxs, file = "cerebxs_QCed.rds")

# -------------------------------------------------------------------------


