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
    
    Obj.list = readRDS(file = "ct_dataset_01.rds")
    Obj.list$min_10_imputate_log2 = Obj.list$filter

    for(i in 1:dim(Obj.list$min_10_imputate_log2)[1]){
        Obj.list$min_10_imputate_log2[i,] = impute_func2(Obj.list$min_10_imputate_log2[i,])
      }
    row.names(Obj.list$min_10_imputate_log2) = row.names(Obj.list$filter)
    saveRDS(Obj.list,file = "ct_dataset_01.rds")
    
}

