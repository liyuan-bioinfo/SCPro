library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)

library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("")

# Imputating missing values with 0.1 * minimum intensity of each protein
{
    impute_func2 <- function(flag.tc_lfq_order){
    min.num <- min(flag.tc_lfq_order[(flag.tc_lfq_order!=0)])/10
    flag.tc_lfq_order[(flag.tc_lfq_order==0)]=min.num
    
    return(log2(flag.tc_lfq_order))
  }

    Obj.list = readRDS(file="ct_dataset_02.rds")
    Obj.list$min_10_impute_log2 = Obj.list$filter
    for(i in 1:dim(Obj.list$min_10_imputate_log2)[1]){#2516
        Obj.list$min_10_impute_log2[i,] = impute_func2(Obj.list$min_10_impute_log2[i,])
    }  
    
}

# Significance calculate
{
    rm(list=ls())
    Obj.list = readRDS("ct_dataset_02.rds")
    
    run_single_comparision = function(groupA,groupB,known_marker){
      temp_meta <- Obj.list$meta
      temp_meta = temp_meta[temp_meta$group %in% c(groupA,groupB),]
      temp_meta$group = factor(temp_meta$group,levels = c(groupA,groupB))
      temp_meta$cell_classification = NULL
      temp_meta = temp_meta[order(temp_meta$group),]
      
      temp.df = Obj.list$no_imputate_log2
      temp.df = temp.df %>% select(temp_meta$sample_id)
      temp.df = temp.df[apply(temp.df,1,sum)!=0,]
      temp.df = Obj.list$min_10_imputate_log2[row.names(temp.df),temp_meta$sample_id]
      
      temp_samples_A = temp_meta$sample_id[which(temp_meta$group==groupA)]
      temp_samples_B = temp_meta$sample_id[which(temp_meta$group==groupB)]
      temp_group_A = groupA
      temp_group_B = groupB
      
      temp_meta.select = temp_meta %>% filter(group==temp_group_A | group==temp_group_B)
      temp_meta.select$group = factor(temp_meta.select$group, levels = c(temp_group_A,temp_group_B))
      temp.df.select = temp.df %>% select(temp_samples_A,temp_samples_B) #df
      
      temp.design <- model.matrix(~temp_meta.select$group)
      colnames(temp.design) <- levels(temp_meta.select$group)
      rownames(temp.design) <- temp_meta.select$sample_id
            
      fit <- lmFit(temp.df.select, temp.design)
      fit <- eBayes(fit, trend=TRUE)
      temp.result.limma <- topTable(fit,  coef=2,n=Inf)
      
      temp.df.select = temp.df.select[row.names(temp.result.limma),]
      temp.df.select$P.Value = temp.result.limma$P.Value
      temp.df.select$fdr = temp.result.limma$adj.P.Val
      temp.df.select$logFC = temp.result.limma$logFC
      temp.df.select$`-log10_P.Value` = -log10(temp.df.select$P.Value)
      temp.df.select$gene = Obj.list$annotation[row.names(temp.df.select),"genename"]
      
      temp.df.select$gene_plot = c()
      temp_ind = temp.df.select$gene %in% known_marker
      temp.df.select$gene_plot[temp_ind] = temp.df.select[temp_ind,"gene"]
      
      return(temp.df.select)
    }
        
    dep_df = run_single_comparision("KPC_425KN","KPC_425KP",c("Foxp3","Cd4","Ikzf2","Klrg1","Kdelr2","Casp1","Map4k1"))
    Obj.list$dep_df = dep_df
    # write.csv(dep_df,file="write/ct_dataset_02_dep.csv")
}
