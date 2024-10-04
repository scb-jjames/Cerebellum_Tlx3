
#Find Tlx3 levels in different cell types of cerebellum at different dev stages
#Data is from : https://doi.org/10.1016/j.cub.2018.07.062
#Analysis by Rahul Jose
#

library(tidyverse)
library(gdata)
library(Seurat)

# -------------------------------------------------------------------------

cerebxcells5 <- readRDS(file = "cerebxcells5_cluster-id.rds")

# For Tlx3 expression across developmental stages
cerebxcells5$cell_type <- Idents(cerebxcells5)

#obtain aggregate expression by developmental stage
tlx3_avgexp2 <-  AggregateExpression(cerebxcells5, features = "Tlx3", group.by = c("ident", "stage"))$RNA
tlx3_avg <- as.data.frame(tlx3_avgexp2)
tlx3 <- tlx3_avg |> 
  pivot_longer(cols = colnames(tlx3_avg), names_to = "cluster", values_to = "count")
tlx3 <- separate(tlx3, col = "cluster", into=c("cell", "stage"), sep = "_")
tlx3 <- filter(tlx3, !cell %in% others)
write.csv(x = tlx3, file = "tlx3_exp.csv")

#For grey background for ggplot
tlx3_grey <- data.frame(cell=as.character(), stage=as.character())
for(x in unique(tlx3$cell)){
  for(y in unique(tlx3$stage)){
    temp <- data.frame(cell=c(x), stage=c(y))
    tlx3_grey <- rbind(tlx3_grey, temp)
  }
}

#Plot
ggplot()+
  geom_tile(data= tlx3_grey, mapping = aes(x=stage, y=cell), fill="gray")+
  geom_tile(data = tlx3, mapping = aes(
    x=factor(stage, levels = c('E10','E11','E12','E13','E14',"E15","E16","E17","P0","P4","P7","P10")),
    y=factor(cell, levels = rev(c("Progenitors", "Granule neuron progenitors", "Granule neurons", "Purkinje cells","Glutamergic CN", 
                                  "GABAergic progenitors", "GABAergic interneurons", "Glial progenitor", "Bergmann glia", "Oligodendrocytes") )), fill=count)
  )+
  scale_fill_continuous(type = "viridis")+
  labs(
    x = "Developmental stages",
    y = "Cell types",
    fill = "Avg. Tlx3 expression",
    title = "Expression of Tlx3 across development in cerebellum"
  ) +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 10, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.title.x = element_text(size = 12, color="black"),
        axis.title.y = element_text(size = 12, color="black"), 
        legend.text = element_text(size = 10, color="black"))
ggsave(filename = "Clusters/tlx3_expr2.jpeg",
       plot = last_plot(),
       dpi = 600)

# -------------------------------------------------------------------------
#Gene expression of other genes across developmental stages

others <- c( "Midbrain cells", "Endothelial cells", "Roof plate cells",  "Microglia", "Meninges", "Erythrocytes", "Ciliated cells")

#function
heatmapx <- function(gene_expr){
  
  #obtain aggregate expression by developmental stage
  gene_expr_avgexp2 <-  AggregateExpression(cerebxcells5, features = gene_expr, group.by = c("ident", "stage"))$RNA
  gene_expr_avg <- as.data.frame(gene_expr_avgexp2)
  gene_expr <- gene_expr_avg |> 
    pivot_longer(cols = colnames(gene_expr_avg), names_to = "cluster", values_to = "count")
  gene_expr <- separate(gene_expr, col = "cluster", into=c("cell", "stage"), sep = "_")
  gene_expr <- filter(gene_expr, !cell %in% others)
  
  #For grey background for ggplot
  gene_expr_grey <- data.frame(cell=as.character(), stage=as.character())
  for(x in unique(gene_expr$cell)){
    for(y in unique(gene_expr$stage)){
      temp <- data.frame(cell=c(x), stage=c(y))
      gene_expr_grey <- rbind(gene_expr_grey, temp)
    }
  }
  
  #plot
  ggplot()+
    geom_tile(data= gene_expr_grey, mapping = aes(x=stage, y=cell), fill="gray")+
    geom_tile(data = gene_expr, mapping = aes(
      x=factor(stage, levels = c('E10','E11','E12','E13','E14',"E15","E16","E17","P0","P4","P7","P10")),
      y=factor(cell, levels = rev(c("Progenitors", "Granule neuron progenitors", "Granule neurons", "Purkinje cells","Glutamergic CN", 
                                    "GABAergic progenitors", "GABAergic interneurons", "Glial progenitor", "Bergmann glia", "Oligodendrocytes") )), fill=count)
    )+
    scale_fill_continuous(type = "viridis")+
    labs(
      x = "Developmental stages",
      y = "Cell types",
      fill = paste("Avg. ", Gene, " expression"),
      title = paste("Expression of ", Gene, " across development in cerebellum")
    ) +
    theme_minimal()+
    theme(axis.text.x = element_text(size = 10, color="black"),
          axis.text.y = element_text(size = 12, color="black"),
          axis.title.x = element_text(size = 12, color="black"),
          axis.title.y = element_text(size = 12, color="black"), 
          legend.text = element_text(size = 10, color="black")
    )
  ggsave(filename = paste("Clusters/", Gene,"_expr2.jpeg"),
         plot = last_plot(), units = "in", width = 8.85, height = 3.47,
         dpi = 600)
  
}

#Plot for other genes
for(x in c("Cadm1",
           "Bbip1",
           "Insc",
           "Cdk2ap1")){
  heatmapx(x)
}

# End ---------------------------------------------------------------------


