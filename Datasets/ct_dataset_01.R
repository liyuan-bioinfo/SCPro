library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)

library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(limma)

setwd("")
# Imputing missing values with 0.1 * minimum intensity of each protein
{
    impute_func2 <- function(flag.tc_lfq_order){
        min.num <- min(flag.tc_lfq_order[(flag.tc_lfq_order!=0)])/10
        flag.tc_lfq_order[(flag.tc_lfq_order==0)]=min.num    
        return(log2(flag.tc_lfq_order))
  }
    
    Obj.list = readRDS(file = "ct_dataset_01.rds")
    Obj.list$min_10_impute_log2 = Obj.list$filter

    for(i in 1:dim(Obj.list$min_10_impute_log2)[1]){
        Obj.list$min_10_impute_log2[i,] = impute_func2(Obj.list$min_10_impute_log2[i,])
      }
    row.names(Obj.list$min_10_impute_log2) = row.names(Obj.list$filter)
    saveRDS(Obj.list,file = "ct_dataset_01.rds")
    
}

# Significiance for 14 cell types
{
    rm(list=ls())
    Obj.list = readRDS(file = "ct_dataset_01.rds")
    
    run_OnevsRest = function(temp_meta, temp_df){
      temp.output.list = list()
      for(i in levels(temp_meta$CellType)[-1]){
        temp_samples_A = temp_meta$SampleID[which(temp_meta$CellType==levels(temp_meta$CellType)[1])]
        temp_samples_B = temp_meta$SampleID[which(temp_meta$CellType==i)]
        temp_group_A = levels(temp_meta$CellType)[1]
        temp_group_B = i
        
        temp_meta.select = temp_meta %>% filter(CellType==temp_group_A | CellType==temp_group_B)
        temp_meta.select$CellType = factor(temp_meta.select$CellType, levels = c(temp_group_A,temp_group_B))
        temp.df.select = temp_df %>% select(temp_samples_A,temp_samples_B) #df
        
        # temp.design <- model.matrix(~0+factor(temp_meta$CellType))
        temp.design <- model.matrix(~temp_meta.select$CellType)
        colnames(temp.design) <- levels(temp_meta.select$CellType)
        rownames(temp.design) <- temp_meta.select$SampleID
        
        
        fit <- lmFit(temp.df.select, temp.design)
        fit <- eBayes(fit, trend=TRUE)
        temp.result.limma <- topTable(fit,  coef=2,n=Inf)
        
        temp.df.select = temp.df.select[row.names(temp.result.limma),]
        temp.df.select$P.Value = temp.result.limma$P.Value
        temp.df.select$fdr = temp.result.limma$adj.P.Val
        temp.df.select$logFC = temp.result.limma$logFC
        
        temp.df.select.sig = temp.df.select %>% filter(P.Value<0.05 & logFC < -1)                
        temp.output.list[[paste(temp_group_A,"-",temp_group_B,sep="")]] = temp.df.select.sig
      }
      ## all-diff gene
      temp.x=row.names(temp.output.list[[1]])
      for(i in (2:length(temp.output.list))){
        temp.x = intersect(temp.x,row.names(temp.output.list[[i]]))
        print(length(temp.x))
      }
      
      if(length(temp.x)==0){
        return(0)
      }
      
      ## average FC
      temp.df.sig = data.frame(pid=temp.x,stringsAsFactors = F)
      temp.logFC = temp.output.list[[1]]
      temp.logFC = temp.logFC[temp.df.sig$pid,"logFC"]
      for(i in (2:length(temp.output.list))){
        t = temp.output.list[[i]][temp.df.sig$pid,"logFC"]
        temp.logFC = temp.logFC+t
      }
      temp.df.sig$logFC = abs(temp.logFC)
      temp.df.sig = temp.df.sig[order(temp.df.sig$logFC,decreasing = T),]
      temp.df.sig$CellType=levels(temp_meta$CellType)[1]
      return(temp.df.sig)
    }
    
    ### PCC vs Rest
    meta <- Obj.list$meta
    
    temp_group = meta %>% dplyr::select(CellType) %>% unique() 
    row.names(temp_group) = temp_group$CellType
    temp_group$CellType = as.character(temp_group$CellType)
    
    temp_OnevsRest_list = list()
    
    ##> PCC
    temp_subgroup = c("PCC",temp_group[temp_group$CellType %in% c("CAF","T4","T8","B","MYE"),])
    temp_OnevsRest_list[["PCC"]] = factor(temp_subgroup, levels = temp_subgroup)
    
    ##>CAF
    temp_subgroup = c("apCAF",temp_group[temp_group$CellType %in% c("myCAF","iCAF"),])
    temp_OnevsRest_list[["apCAF"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("myCAF",temp_group[temp_group$CellType %in% c("apCAF","iCAF"),])
    temp_OnevsRest_list[["myCAF"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("iCAF",temp_group[temp_group$CellType %in% c("apCAF","myCAF"),])
    temp_OnevsRest_list[["iCAF"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("CAF",temp_group[temp_group$CellType %in% c("PCC","T4","T8","B","MYE"),])
    temp_OnevsRest_list[["CAF"]] = factor(temp_subgroup, levels = temp_subgroup)
    
    ##>Lymphoid
    temp_subgroup = c("T4",temp_group[temp_group$CellType %in% c("T8","B"),])
    temp_OnevsRest_list[["T4"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("T8",temp_group[temp_group$CellType %in% c("T4","B","Treg"),])
    temp_OnevsRest_list[["T8"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("Treg",temp_group[temp_group$CellType %in% c("T4","T8","B"),])
    temp_OnevsRest_list[["Treg"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("B",temp_group[temp_group$CellType %in% c("T4","T8","Treg"),])
    temp_OnevsRest_list[["B"]] = factor(temp_subgroup, levels = temp_subgroup)
    
    ##> MYE
    temp_subgroup = c("MYE",temp_group[temp_group$CellType %in% c("PCC","T4","T8","B"),])
    temp_OnevsRest_list[["MYE"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("NEU",temp_group[temp_group$CellType %in% c("MO","MAC","DC"),])
    temp_OnevsRest_list[["NEU"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("MO",temp_group[temp_group$CellType %in% c("NEU","MAC","DC"),])
    temp_OnevsRest_list[["MO"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("MAC",temp_group[temp_group$CellType %in% c("MO","NEU","DC"),])
    temp_OnevsRest_list[["MAC"]] = factor(temp_subgroup, levels = temp_subgroup)
    temp_subgroup = c("DC",temp_group[temp_group$CellType %in% c("MO","MAC","NEU"),])
    temp_OnevsRest_list[["DC"]] = factor(temp_subgroup, levels = temp_subgroup)    
    
    temp_out.df = data.frame()
    for(i in temp_OnevsRest_list){
      temp_submeta = meta %>% filter(CellType %in% i)
      temp_submeta$CellType = factor(x=temp_submeta$CellType, levels = i)
      temp_submeta = temp_submeta[order(temp_submeta$CellType),]
      
      temp_subdf = Obj.list$min_10_impute_log2 %>% 
        select(temp_submeta$SampleID)
      temp.out = run_OnevsRest(temp_submeta, temp_subdf)
      if(temp.out!=0){
        temp_out.df = rbind(temp_out.df,temp.out)  
      }
    }
    
    temp_out.df$genename = Obj.list$annotation[temp_out.df$pid,"genename"]
    temp_out.df[temp_out.df$pid %in% Obj.list$uniprot_mouse_loci$Entry,"Loca"]="Surface"
            
    Obj.list$dep_df = temp_out.df
    saveRDS(Obj.list,file = "ct_dataset_01.rds")
    write.csv(temp_out.df,file="write/OnevsRest_DEP_Min10Impute.csv"),row.names = F,na = "")
  }
