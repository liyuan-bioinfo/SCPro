library(readxl)
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggplot2)

setwd("")

# for Figure-2J, cluster
{
  data.df = readxl::read_xlsx("combined_protein_20230822.xlsx",sheet=1) %>% as.data.frame() # output from Perseus with one-way ANOVA
  names(data.df)=gsub(names(data.df),pattern = " MaxLFQ.*$",replacement = "")
  data.df = data.df[which(data.df$`ANOVA Significant` == "+"),]#1776
  
  meta.df = data.frame(SampleID=names(data.df)[3:14])
  meta.df$group = gsub(meta.df$SampleID,pattern = "[1-4]",replacement = "")  
  row.names(data.df) = data.df$`Protein ID`
  
  p1=pheatmap::pheatmap(data.df[,meta.df$SampleID],show_rownames = T,scale="row",fontsize = 2,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  pdf("write/pheatmap_cluster_showpid.pdf",width=6,height = 60)
  print(p1)
  dev.off()
  
  temp.df = data.frame(pid=p1$tree_row$labels[p1$tree_row$order])
  write.csv(temp.df,file="write/pheatmap_cluster_table.csv")
}

# for figure-2J, GO-BP annotation
{  
  pheatmap_cluster.df = read.csv(file="write/pheatmap_cluster_table.csv",header = T,row.names = 1)
  pheatmap_cluster.df$group = factor(pheatmap_cluster.df$group, levels = c("Acinar", "Tumor", "Lymph"))
  
  tran.df=clusterProfiler::bitr(pheatmap_cluster.df$pid,fromType = "UNIPROT",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
  pheatmap_cluster.df = pheatmap_cluster.df %>% merge(tran.df,by.x="pid",by.y="UNIPROT",all.x=T)    
  formula_res <- compareCluster(data = pheatmap_cluster.df,ENTREZID~group,fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,
                                ont = "ALL", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE                                
  )
  
  write.csv(summary(formula_res), file="write/figure_2J_GOBP_enrichment.csv")  
}
