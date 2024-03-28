library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)

library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("")
# Signficance Analysis with one-way ANOVA, for Figure 4c, 4d and 4e.
{
    rm(list=ls())
    Obj.list = readRDS(file="sp_dataset_03.rdata")
    regions = c("Acinar","PanIN","PDAC")
    meta.df = Obj.list$meta %>% dplyr::filter(Group %in% regions)
        
    meta.df$Regions = factor(meta.df$Group,levels = regions)
    meta.df = meta.df %>% arrange(Regions)
    input.df = Obj.list$impute[,meta.df$SampleID]
    
    Pvalue = c()
    for(i in 1:dim(input.df)[1]){
        pid.df = input.df[i,] %>% t() %>% as.data.frame()
        names(pid.df) = "pid"
        pid.df$SampleId = row.names(pid.df)
        pid.df$Region = meta.df$Regions
        p.value = summary(aov(data=pid.df,pid~Region))[[1]]$`Pr(>F)`[1]
        Pvalue = c(Pvalue,p.value)
        
    }
    input.df$Pvalue = Pvalue
    input.df$fdr = p.adjust(Pvalue,method = "BH")
    
    median.df = Obj.list$impute_median
    median.df = median.df %>% dplyr::select(Acinar,PanIN,PDAC)
    
    dep.df = data.frame()
    for(i in 1:3){
        temp.df = median.df[,i] - median.df[,1:3]
        temp_pid = names(which(rowSums(temp.df > log2(2)) == 2)) # The fold change with one vs the rest strategy. For figure 4e, should change to log2(1.2)
        temp.df2 = data.frame(pid=temp_pid)
        temp.df2$region = names(median.df)[i]
        temp.df2$log2FC = apply(temp.df[temp_pid,-i],1,median)
        dep.df = rbind(dep.df,temp.df2)
        
    }
    row.names(dep.df) = dep.df$pid
    dep.df$genename = Obj.list$anno[dep.df$pid,"Gene"]
    
    temp.df = input.df %>% dplyr::select(Pvalue, fdr)
    temp.df$pid = row.names(temp.df)
    dep.df = dep.df %>% merge(temp.df,by="pid",all.x=T,all.y=F)#add sig
    
    dep.df$region = factor(dep.df$region,levels=regions)
    dep.df = dep.df %>% arrange(region) # Table 1
    
    Obj.list$PDAC_stage = dep.df
    saveRDS(Obj.list,file = "sp_dataset_03.rdata")    
}

# Signficance Analysis with two-sided Student's t-test, for Figure 4g and 4h.
{
    rm(list=ls())
    Obj.list = readRDS(file="write/sp_dataset_03.rdata")
    regions = c("LN","IT")
    meta.df = Obj.list$meta %>% dplyr::filter(Group %in% regions)
        
    meta.df$Regions = factor(meta.df$Group,levels = regions)
    meta.df = meta.df %>% arrange(Regions)
    input.df = Obj.list$filter[,meta.df$SampleID]
    input.df = input.df[-which(rowSums(input.df)==0),]
    
    Pvalue = c()
    log2FC = c()
    for(i in 1:dim(input.df)[1]){
        pid.df = input.df[i,] %>% t() %>% as.data.frame()
        names(pid.df) = "pid"
        pid.df$SampleId = row.names(pid.df)
        pid.df$Region = meta.df$Regions
        p.value = t.test(data=pid.df,pid~Region)$p.value
        Pvalue = c(Pvalue,p.value)
        
        temp_FC = median(pid.df$pid[1:3]) - median(pid.df$pid[4:6])
        log2FC = c(log2FC, temp_FC)
        
    }
    input.df$Pvalue = Pvalue
    input.df$fdr = p.adjust(Pvalue,method = "BH")
    input.df$log2FC = log2FC
    
    input.df$pid = Obj.list$anno[row.names(input.df),"pid"]
    input.df$genename = Obj.list$anno[row.names(input.df),"Gene"]
    
    Obj.list$Lym_DEP = input.df
    saveRDS(Obj.list,file = "sp_dataset_03.rdata")
}
