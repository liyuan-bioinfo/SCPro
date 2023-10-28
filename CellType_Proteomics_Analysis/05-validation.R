#'@author LiYuan
#'@desc Cd25 vs Klrg1 另外的样本、重新搜库
#'@time 20221010
#'@desc Fig4 
#'

rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(ggbiplot)
library(corrplot)
library(limma)
library(RColorBrewer)
library(ggrepel)
library(org.Mm.eg.db)
setwd("E:\\XYF（许燕芬）\\FACS_SISPROT\\workspace\\20221128_final\\05-validation/")
options(stringsAsFactors=F)

### pre-treat, filter, log2, imputate
{    
  ## MSFragger
  Obj.list = list()
  Obj.list$raw = read.csv(file="KLRG1_allsubtypes_combined_protein(1).tsv",stringsAsFactors = F,
                                      fill = FALSE,sep="\t")
  
  Obj.list$raw = Obj.list$raw %>% filter(Organism=="Mus musculus OX=10090")#5102
  
  ## surface protein 
  # temp.list = readRDS("MMU.pm.protein.database_4518_20220626.rda")
  # Obj.list$uniprot_mouse_loci = temp.list$combined_pm.db
  
  
  sample_ids = names(Obj.list$raw[1,c(208:231)])
  
  ##> annotation table2: id, group
  temp.meta = data.frame(raw_sample_id=sample_ids)
  temp.meta$sample_id = gsub(x=temp.meta$raw_sample_id,replacement = "",pattern = "_[0-9]+.MaxLFQ.Total.Intensity")
  temp.meta$group = gsub(temp.meta$sample_id,replacement = "",pattern = "_[1-3]$")
  temp.meta$group = factor(temp.meta$group,level=c("KPC_425KN","KPC_425KP",
                                                   "KPC_425N","KPC_425P",
                                                   "KPC_4KN","KPC_4KP",
                                                   "KPC_8KN","KPC_8KP"))
  row.names(temp.meta) = temp.meta$sample_id
  temp.meta = temp.meta[order(temp.meta$group),]
  
  temp.meta$group_color = factor(
                                   c(rep("#08306B",3),rep("#08519C",3),rep("#2171B5",3),rep("#4292C6",3),
                                   rep("#FD8D3C",3),rep("#F16913",3),
                                   rep("#ADDD8E",3),rep("#78C679",3)),
                                 levels=c( 
                                            "#08306B","#08519C","#2171B5","#4292C6",#CAF
                                            "#FD8D3C","#F16913",
                                            "#ADDD8E","#78C679"))
  Obj.list$meta = temp.meta
  # temp.meta$group_color = factor(c(rep("#6BAED6",3),rep("#004529",3)),
  #                                levels=c(  "#6BAED6","#004529"))
  write.csv(file="test.csv",df)
  
  
  Obj.list$filter = Obj.list$raw[,temp.meta$raw_sample_id] 
  names(Obj.list$filter) = temp.meta$sample_id
  row.names(Obj.list$filter) = Obj.list$raw$Protein.ID
  
  Obj.list$quantified = Obj.list$raw[,temp.meta$raw_sample_id] 
  names(Obj.list$quantified) = temp.meta$sample_id
  row.names(Obj.list$quantified) = Obj.list$raw$Protein.ID
  Obj.list$quantified = Obj.list$quantified[apply(Obj.list$quantified,1,sum)!=0,]
  
  Obj.list$identified = Obj.list$raw[,c(136:159)] 
  row.names(Obj.list$identified) = Obj.list$raw$Protein.ID
  
  Obj.list$identified = Obj.list$identified[apply(Obj.list$identified,1,sum)!=0,]
  names(Obj.list$identified) = gsub(x=names(Obj.list$identified),replacement = "",pattern = "_[0-9]+.Total.Intensity")
  Obj.list$identified = Obj.list$identified[,Obj.list$meta$sample_id]
  
  
  ##> annotation table1:pid, genename
  Obj.list$annotation = data.frame(pid=Obj.list$raw$Protein.ID,
                                                   genename=Obj.list$raw$Gene)
  row.names(Obj.list$annotation) = Obj.list$annotation$pid


  ##> at least two valid values in one group
  # Obj.list$filter = Obj.list$filter[which(apply(Obj.list$filter, 1,sum)!=0),]
  temp.valid = c()
  for(i in 1:dim(Obj.list$filter)[1]){
    temp.df = Obj.list$filter[i,]
    for(j in levels(temp.meta$group)){
      temp.id = temp.meta[temp.meta$group%in% j,"sample_id"]
      if(length(which(temp.df[,temp.id]>0))>=2){
        temp.valid = c(temp.valid,i)
        break
      }
    }
  }
  Obj.list$filter = Obj.list$filter[temp.valid,]
  
  ## imputate 1: no imputate
  Obj.list$no_imputate_log2 = log2(Obj.list$filter+1)
  row.names(Obj.list$no_imputate_log2) = row.names(Obj.list$filter)
  
  impute_func1 <- function(flag.tc_lfq_order, shift = 3, width = 0.3)
  {
    NA.num <- sum(is.na(flag.tc_lfq_order))
    vec_RemoveNA <- as.vector(as.matrix(flag.tc_lfq_order))
    vec_RemoveNA <- na.omit(vec_RemoveNA)
    
    missingValueMean <- mean(vec_RemoveNA,na.rm = TRUE) - shift*sd(vec_RemoveNA,na.rm = TRUE)
    missingVlaueSd <- width*sd(vec_RemoveNA,na.rm = TRUE)
    set.seed(1234)
    NA.replace.num <- rnorm(NA.num, mean = missingValueMean, sd=missingVlaueSd)
    flag.tc_lfq_order[is.na(flag.tc_lfq_order)] <- NA.replace.num
    return(flag.tc_lfq_order)
  }
  
  ##> imputate2: random distribution 
  Obj.list$mann_imputate_log2 = log2(Obj.list$filter)
  for(i in 1:dim(Obj.list$mann_imputate_log2)[1]){
    Obj.list$mann_imputate_log2[i,] = as.numeric(gsub(pattern = "-Inf",NA,x=Obj.list$mann_imputate_log2[i,]))
    Obj.list$mann_imputate_log2[i,] = impute_func1(Obj.list$mann_imputate_log2[i,])
  }
  row.names(Obj.list$mann_imputate_log2) = row.names(Obj.list$filter)
  # write.csv(Obj.list$mann_imputate_log2,"test.csv")
  ###> imputate3: min/10
  impute_func2 <- function(flag.tc_lfq_order)
  {
    min.num <- min(flag.tc_lfq_order[(flag.tc_lfq_order!=0)])/10
    flag.tc_lfq_order[(flag.tc_lfq_order==0)]=min.num
    
    return(log2(flag.tc_lfq_order))
  }
  Obj.list$min_10_imputate_log2 = Obj.list$filter
  for(i in 1:dim(Obj.list$min_10_imputate_log2)[1]){#2516
    Obj.list$min_10_imputate_log2[i,] = impute_func2(Obj.list$min_10_imputate_log2[i,])
  }
  row.names(Obj.list$min_10_imputate_log2) = row.names(Obj.list$filter)
  
  #row_means
  temp.samples.anno = Obj.list$meta 
  temp.mean.df = data.frame(row.names = row.names(Obj.list$filter))
  for(i in levels(temp.samples.anno$group)){
    temp_index = which(temp.samples.anno$group==i)  
    temp_df = Obj.list$min_10_imputate_log2[,temp_index]
    temp.mean.df[,i] = rowMeans(temp_df,na.rm = TRUE)
  }
  #zscore
  temp.mean.zscore=apply(temp.mean.df,1,function(x){ 
    (x-mean(x))/sd(x)
  })
  temp.mean.zscore = temp.mean.zscore %>% t() %>% as.data.frame()
  
  Obj.list$min_10_imputate_log2_row_means = temp.mean.df
  Obj.list$min_10_imputate_log2_row_means_zscore = temp.mean.zscore
  
  #row_means
  temp.samples.anno = Obj.list$meta 
  temp.mean.df = data.frame(row.names = row.names(Obj.list$filter))
  for(i in levels(temp.samples.anno$group)){
    temp_index = which(temp.samples.anno$group==i)  
    temp_df = Obj.list$mann_imputate_log2[,temp_index]
    temp.mean.df[,i] = rowMeans(temp_df,na.rm = TRUE)
  }
  #zscore
  temp.mean.zscore=apply(temp.mean.df,1,function(x){ 
    (x-mean(x))/sd(x)
  })
  temp.mean.zscore = temp.mean.zscore %>% t() %>% as.data.frame()
  
  Obj.list$mann_imputate_log2_row_means = temp.mean.df
  Obj.list$mann_imputate_log2_row_means_zscore = temp.mean.zscore
  
  saveRDS(Obj.list,file = paste0("Validation_Cd25_Klrg1_All_Obj_",Sys.Date(),".rds"))
  
}
## save table
{
  Obj.list = readRDS(file="Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")
  temp.df = Obj.list$quantified
  temp.df$genename = Obj.list$annotation[row.names(temp.df),"genename"]
  write.csv(x = temp.df,
            file=paste("Quantified_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
  temp.df = Obj.list$identified
  temp.df$genename = Obj.list$annotation[row.names(temp.df),"genename"]
  write.csv(x = temp.df,
            file=paste("Identified_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
  
}

#------------------------------------------------------------------------
## Plot
{
  
  #Fig2,C 
  ## Bar plot showing the protein numbers of quantified and identified proteins 
  {
    rm(list=ls())
    Obj.list = readRDS(file = "../05-validation/Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")
    temp.df = Obj.list$identified %>% t() %>% as.data.frame()
    
    meta.df = Obj.list$meta
    meta.df$SampleID = meta.df$sample_id
    meta.df$CellType = meta.df$group
    meta.df = meta.df %>% select(SampleID,CellType) %>% unique()
    # row.names(plot.data) = plot.data$group
    plot.data = meta.df %>% select(CellType) %>% unique()
    row.names(plot.data) = plot.data$CellType
    
    temp_sum = c()
    for(i in levels(plot.data$CellType)){
      temp_index = which(meta.df$CellType==i)  
      tt = temp.df[temp_index,]
      ttt = apply(tt,2,function(x){
        return(sum(x,na.rm = TRUE))
      })
      temp_len = length(which(ttt >0))
      temp_sum = c(temp_sum,temp_len)
      
    }
    plot.data$protein_sum =temp_sum
    
    ## quantified
    plot.data2 = meta.df %>% select(CellType) %>% unique()
    row.names(plot.data2) = plot.data2$CellType
    temp.df2 = Obj.list$quantified %>% t() %>% as.data.frame()
    
    temp_sum = c()
    for(i in levels(plot.data2$CellType)){
      temp_index = which(meta.df$CellType==i)  
      # t = unique(which(temp.df[temp_index,]!=0))
      tt = temp.df2[temp_index,]
      ttt = apply(tt,2,function(x){
        return(sum(x,na.rm = TRUE))
      })
      temp_len = length(which(ttt >0))
      temp_sum = c(temp_sum,temp_len)
      
    }
    plot.data2$protein_sum =temp_sum
    
    #combine
    plot.data["Total","protein_sum"]=dim(Obj.list$identified)[1]
    plot.data$CellType = row.names(plot.data)
    plot.data$protein_loca = "Identified"
    
    plot.data2["Total","protein_sum"]=dim(Obj.list$quantified)[1]
    plot.data2$CellType = row.names(plot.data2)
    plot.data2$protein_loca = "Quantified"
    
    plot.data3 = plot.data
    plot.data3$protein_sum = plot.data$protein_sum - plot.data2$protein_sum
    plot.data.bar = rbind(plot.data3,plot.data2)
    plot.data.bar$CellType = factor(x = plot.data.bar$CellType,levels = c(levels(meta.df$CellType),"Total"))
    
    temp.p1 <- ggplot(plot.data.bar, aes(x=CellType,y=protein_sum,fill=protein_loca)) + 
      geom_bar(stat = "identity",width=0.8)+
      scale_fill_grey(start = 0.75, end = 0.5) +
      theme_classic()+
      theme(        
        plot.title = element_text(size=16,hjust = 0.5),
        axis.title.y.left  = element_text(size=16),
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1,color = c(levels(Obj.list$meta$group_color),"Black")),
        axis.ticks.x = element_blank(),axis.line.x.bottom = element_blank(),
        legend.position=c(0.5,0.8),legend.direction = "horizontal",legend.text = element_text(size=12),
        legend.background = element_blank(),panel.grid = element_blank(),
        plot.background = element_blank(),strip.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank()
      )+
      labs(y="Protein groups(number)",fill="",x="")+
      # ylab("Protein groups(number)")+
      scale_y_continuous(expand = c(0,0),limits=c(0,6000),breaks=c(0,2000,4000,5000,6000),label=c(0,2000,4000,5000,6000))
    
    
    ## Output  
    plot.data.bar2 = rbind(plot.data,plot.data2)
    write.csv(plot.data.bar2,
              file = paste0("Fig2_C_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".csv"),
              fileEncoding = "UTF-8",row.names = F,quote = TRUE)
    
    pdf(paste0("Fig2_C_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".pdf"),
        width=6,height = 4)
    print(temp.p1)
    dev.off()
  }   
  #box plot of each samples
  {
    {
      rm(list=ls())
      Obj.list = readRDS(file = "Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")
      # all proteins
      temp.df = Obj.list$identified %>% t() %>% as.data.frame()
      temp = apply(temp.df,MARGIN = 1,function(x){
        return(length(which(x!=0)))
      })
      plot.data = Obj.list$meta
      plot.data$protein_num = temp
      plot.data$protein_loca = "Identified"
      
      temp.df2 = Obj.list$quantified %>% t() %>% as.data.frame()
      temp2 = apply(temp.df2,MARGIN = 1,function(x){
        return(length(which(x!=0)))
      })
      plot.data2 = plot.data
      plot.data2$protein_num = temp2
      plot.data2$protein_loca = "Quantified"
      
      plot.data.box = rbind(plot.data, plot.data2)
      
      temp.p1 = plot.data.box %>%
        ggplot( aes(x=group, y=protein_num, fill=protein_loca)) +
        geom_boxplot() +
        geom_jitter(color="black", size=1, alpha=10) +
        # stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="red", fill="red") +
        theme_bw()+
        theme(
          #legend.position=c("bottom"),legend.box = "horizontal",
          plot.title = element_text(size=16,hjust = 0.5),
          axis.text = element_text(size=12),
          axis.title = element_text(size=18),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     color = levels(Obj.list$meta$group_color)),
          legend.background = element_blank(),panel.grid = element_blank(),
          plot.background = element_blank(),strip.background = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank()
        ) + labs(fill="")+
        scale_fill_grey(start = 0.75, end = 0.25) +
        labs(title="All identified and quantified proteins",x="Cell types", y = "Protein number")
      
      ## save data
      write.csv(plot.data.box,
                file = paste0("FigS2_C_ProteinNumber_PG_IdentifiedvsQuantified_Boxplot_",Sys.Date(),".csv"))
      # print(temp.p1)    
      pdf(paste0("FigS2_C_ProteinNumber_PG_IdentifiedvsQuantified_Boxplot_",
                 Sys.Date(),".pdf"), width=6,height = 4)
      print(temp.p1)
      dev.off()
      
    }
    
  }
#expression distribution
  {
    rm(list=ls())
    Obj.list = readRDS(file = "Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")
    
    temp.df = Obj.list$mann_imputate_log2 %>%
      gather(key="name",value="value")
    meta = Obj.list$meta
    meta$SampleID = meta$sample_id
    meta$CellType = meta$group
    
    temp.df = temp.df %>% merge(meta,by.x="name",by.y="SampleID")
    temp.df$name = factor(temp.df$name,levels=meta$SampleID)
    
    temp.p1=ggplot(data = temp.df,mapping = aes(x=name,y=value,fill=group_color)) + 
      geom_boxplot(fill=Obj.list$meta$group_color,outlier.colour = "black",color="black",outlier.size = 0.1) +
      # geom_jitter(color="black", size=1, alpha=1) +
      theme_bw()+
      theme(
        legend.position="none",
        plot.title = element_text(size=18,hjust = 0.5),axis.text = element_text(size=10),
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank()
      ) +
      stat_boxplot(geom = "errorbar", width = 0.15) +
      labs(title=paste("Expression distribution of all samples"),x="Sample ID", y = "Log2(LFQ.Intensity")
    
    
    pdf(file = paste0("ExpressionDistribution_MannImputated_Boxplot_",Sys.Date(),".pdf"),width=6,height=4)
    print(temp.p1)
    dev.off()
    
    write.csv(temp.df,file=paste0("ExpressionDistribution_MannImputated_Boxplot_",Sys.Date(),".csv",sep=""))
  }
#Fig6_C
### 03-PCA
  {
  rm(list=ls())
  Obj.list = readRDS(file = "Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")
  temp.df = Obj.list$min_10_imputate_log2 %>% t() %>% as.data.frame()
  
  temp.df.anno = Obj.list$meta
  temp.df.anno = temp.df.anno[row.names(temp.df),]
  temp.df.anno$sample_id = c(
    paste0("CD4+CD25+Klrg1-_",1:3),paste0("CD4+CD25+Klrg1+_",1:3),
    paste0("CD4+CD25-_",1:3),paste0("CD4+CD25+_",1:3),
    paste0("CD4+Klrg1-_",1:3),paste0("CD4+Klrg1+_",1:3),
    paste0("CD8+Klrg1-_",1:3),paste0("CD8+Klrg1+_",1:3)
  )
  temp.df.anno$group =  gsub(temp.df.anno$sample_id,replacement = "",pattern = "_.*$")
  temp.df.anno$group = factor(temp.df.anno$group,levels = (
    c("CD4+CD25+Klrg1-","CD4+CD25+Klrg1+",
    "CD4+CD25-","CD4+CD25+",
    "CD4+Klrg1-","CD4+Klrg1+",
    "CD8+Klrg1-","CD8+Klrg1+")
    
  ))
  temp.df.anno$shape = factor(      
    c(rep("CD4+CD25+Klrg1",6),
      rep("CD4+CD25",6),
      rep("CD4+Klrg1",6),
      rep("CD8+Klrg1",6)),
    
    levels = (
    c("CD4+CD25+Klrg1",
      "CD4+CD25",
      "CD4+Klrg1",
      "CD8+Klrg1")
    
  ))
  pca_result <- prcomp(temp.df,scale=T,center = T)
  
  df1 = as.data.frame(pca_result$x)
  summ1 <- summary(pca_result)
  xlab1 <- paste0("PC1 (",round(summ1$importance[2,1]*100,2),"%)")
  ylab1 <- paste0("PC2 (",round(summ1$importance[2,2]*100,2),"%)")
  #
  p1=ggplot(data = df1,aes(x = PC1,y = PC2,color = temp.df.anno$group,
                           shape=temp.df.anno$shape
                           ))+
    # stat_ellipse(aes(fill = temp.df.anno$classfication),
    #              type = "norm",geom = "polygon",alpha = 0.25,
    #              color = NA
    # )+ # 添加置信椭圆
    geom_point(size = 3.5)+
    # geom_text_repel(data=df1,aes(x=PC1,y=PC2,color=temp.df.anno$group,label=temp.df.anno$sample_id))+
    labs(x = xlab1,y = ylab1,color = "",title = "",shape="")+
    guides(fill = "none")+
    theme_bw()+
    # scale_fill_manual(values = levels(temp.df.anno$))+
    scale_colour_manual(values = levels(temp.df.anno$group_color))+
    theme(plot.background = element_blank(),legend.background = element_blank(),
          panel.background = element_blank(),panel.grid = element_blank(),
          axis.text = element_text(size = 12),axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
  
  
  pdf(file=paste("1_ProteinExpression_Cd25vsKlrg1_All/Fig6_C_PCA_MinImputated.",Sys.Date(),".pdf",sep=""),
      width = 8,height=5)
  print(p1)
  dev.off()
  write.csv("pca.csv",x = df1)
}

#Fig6_D E
### 07-Comparision analysis,Volcano plot
  {
    rm(list=ls())
    Obj.list = readRDS("Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")
    
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
      
      # temp.design <- model.matrix(~0+factor(temp_meta$group))
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
      # known_marker = c("Cd8a","Cd8b","Cd4")
      temp_ind = temp.df.select$gene %in% known_marker
      temp.df.select$gene_plot[temp_ind] = temp.df.select[temp_ind,"gene"]
      
      logFC_cutoff <- log2(2)
      log10_P_Value_cutoff <- -log10(0.05)
      
      plot1 <- ggplot(data = temp.df.select,aes(x = logFC,y = `-log10_P.Value`))+
        geom_point(data = subset(temp.df.select,abs(logFC)<logFC_cutoff),size=2,
                   col = 'gray',alpha = 0.4)+
        geom_point(data = subset(temp.df.select,abs(`-log10_P.Value`)<log10_P_Value_cutoff & abs(logFC)>logFC_cutoff),
                   size=2,col = 'gray',alpha = 0.4)+
        geom_point(data = subset(temp.df.select,abs(`-log10_P.Value`)>log10_P_Value_cutoff & logFC>logFC_cutoff),
                   size=2,col = 'red',alpha = 1)+
        geom_point(data = subset(temp.df.select,abs(`-log10_P.Value`)>log10_P_Value_cutoff & logFC< -logFC_cutoff),
                   size=2,col = 'blue',alpha = 1)+
        geom_point(data = subset(temp.df.select,!is.na(gene_plot)),
                   size=2,colour="black",alpha = 0.4)+
        
        theme_bw()+
        labs(x='log2(fold change)',y='-log10(p-value)')+
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              plot.title = element_text(hjust=0.5,size = 14),
              text = element_text(size=14),
              legend.position = 'none',
              axis.line = element_line(colour = "black"))+
        
        geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),lty = 3,col = 'black',lwd = 0.4)+
        geom_hline(yintercept = log10_P_Value_cutoff,lty = 3,col = 'black',lwd = 0.4)  +
        # geom_text_repel(data = subset(temp.df.select,abs(`-log10_P.Value`)>log10_P_Value_cutoff & abs(logFC)>logFC_cutoff),
        geom_text_repel(data = temp.df.select,min.segment.length=0,
                        aes(label = gene_plot),size = 7,col = 'black') +
        ggtitle(paste("Volcano Plot of ",groupA," vs ",groupB,sep=""))
      return(list(plot1=plot1,table1=temp.df.select))
    }
    
    output_name = "."
    ##> main figure
    ## Cd4+Cd25+Klrg1- vs Cd4+Cd25+Klrg1+
    comparision_name = "Fig6_E_Cd4+Cd25+Klrg1-vsCd4+Cd25+Klrg1+"
    
    sig_output.obj = run_single_comparision("KPC_425KN","KPC_425KP",c("Foxp3","Cd4","Ikzf2","Klrg1","Kdelr2","Casp1","Map4k1"))
    pdf(file=paste0(output_name,"/",comparision_name,"_VolcanoPlot_",Sys.Date(),".pdf"))
    print(sig_output.obj$plot1)
    dev.off()
    write.csv(sig_output.obj$table1,file=paste0(output_name,"/",comparision_name,"_",Sys.Date(),".csv"))
    
    ##> GSEA
    temp_meta = Obj.list$meta
    temp_meta = temp_meta[temp_meta$group %in% c("KPC_425KN","KPC_425KP"),]
    target_pid = Obj.list$no_imputate_log2[,temp_meta$sample_id]
    target_pid = target_pid[apply(target_pid,1,sum)!=0,]
    
    df = sig_output.obj$table1
    target.df = df[row.names(target_pid),] #%>% filter(logFC>=1 | logFC<=-1) %>% filter(P.Value<=0.05)
    target.df$rank = target.df$`-log10_P.Value` * target.df$logFC
    target.df = target.df[order(target.df$rank,decreasing = T),]
    
    geneList = target.df$rank
    names(geneList) = row.names(target.df)
    re = clusterProfiler::gseGO(geneList =geneList,ont="BP",OrgDb=org.Mm.eg.db,
                                  keyType = "UNIPROT",
                                  by = "fgsea",pvalueCutoff = 1,seed=123
    )
    write.csv(summary(re),file = "GSEA_Klrg1_DEP_GSEA_GO.csv")
    
    enrichplot::gseaplot2(re,geneSetID="GO:0002275")
    enrichplot::gseaplot2(re,geneSetID="GO:0050870")
    clusterProfiler::gseaplot(re,geneSetID="GO:0050870")
    
    gmt = clusterProfiler::read.gmt("m5.go.bp.v2022.1.Mm.entrez.gmt")
    gmt.tran <- bitr(gmt$gene, fromType = "ENTREZID", toType = c("UNIPROT"),OrgDb = org.Mm.eg.db)
    gmt.tran.df = gmt.tran %>% merge(gmt,by.x="ENTREZID",by.y="gene",all.x=FALSE,all.y=FALSE) ##883
    
    gmt.tran_pid_term.df = gmt.tran.df %>% select(UNIPROT,term) %>% unique()
    names(gmt.tran_pid_term.df) = c("gene","term")
    gmt.tran_pid_term.df = gmt.tran_pid_term.df[,2:1]
    re = clusterProfiler::GSEA(geneList =geneList,TERM2GENE = gmt.tran_pid_term.df,
                               # keyType = "UNIPROT",
                               by = "fgsea",pvalueCutoff = 1,seed=123
    )
    
    p1=enrichplot::gseaplot2(re,geneSetID="GOBP_MYELOID_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")
    p2=enrichplot::gseaplot2(re,geneSetID="GOBP_POSITIVE_REGULATION_OF_LYMPHOCYTE_ACTIVATION")
    
    write.csv(summary(re),file = "GSEA_Klrg1_DEP_GSEA.csv")
    pdf("Immune.pdf")
    print(p1)
    dev.off()
    
    pdf("T4.pdf")
    print(p2)
    dev.off()
    ##> figure
    ## Cd4+Cd25+ vs Cd4+Klrg1+
    comparision_name = "FigS6_B_Cd4+Cd25+vsCd4+Klrg1+"
    
    sig_output.obj = run_single_comparision("KPC_425P","KPC_4KP",c("Foxp3","Cd4","Ikzf2","Klrg1"))
    pdf(file=paste0(output_name,"/",comparision_name,"_VolcanoPlot_",Sys.Date(),".pdf"))
    print(sig_output.obj$plot1)
    dev.off()
    write.csv(sig_output.obj$table1,file=paste0(output_name,"/",comparision_name,"_",Sys.Date(),".csv"))
    
    {
    ## Cd4+Klrg1- vs Cd4+Klrg1+
    comparision_name = "Cd4+Klrg1-vsCd4+Klrg1+"
    
    sig_output.obj = run_single_comparision("KPC_4KN","KPC_4KP",c("Foxp3","Cd4","Klrg1","Ikzf2"))
    pdf(file=paste0(output_name,"/",comparision_name,"_VolcanoPlot_",Sys.Date(),".pdf"))
    print(sig_output.obj$plot1)
    dev.off()
    write.csv(sig_output.obj$table1,file=paste0(output_name,"/",comparision_name,"_",Sys.Date(),".csv"))
    
    ## Cd8+Klrg1- vs Cd8+Klrg1+
    comparision_name = "Cd8+Klrg1-vsCd8+Klrg1+"
    
    sig_output.obj = run_single_comparision("KPC_8KN","KPC_8KP",c("Foxp3","Cd4","Klrg1","Cd8a","Cd8b","Cd3"))
    pdf(file=paste0(output_name,"/",comparision_name,"_VolcanoPlot_",Sys.Date(),".pdf"))
    print(sig_output.obj$plot1)
    dev.off()
    write.csv(sig_output.obj$table1,file=paste0(output_name,"/",comparision_name,"_",Sys.Date(),".csv"))
    }
  }

#Fig6_F
  ### 08-GO annotation
  {
    library(DOSE)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(topGO)
    library(clusterProfiler)
    rm(list=ls())
    run_GO_BP_BG = function(sig.pid,bg.pid,items=10){
      sig.pid.tran <- bitr(sig.pid, fromType = "UNIPROT", toType = c("ENSEMBL", "ENTREZID","SYMBOL"), 
                           OrgDb = org.Mm.eg.db)
      bg.pid.tran <- bitr(bg.pid, fromType = "UNIPROT", toType = c("ENSEMBL", "ENTREZID","SYMBOL"), 
                          OrgDb = org.Mm.eg.db)
      
      ego_BP <- enrichGO(gene = sig.pid.tran$ENTREZID, 
                         OrgDb = org.Mm.eg.db,ont = "BP", pAdjustMethod = "BH",
                         pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE,
                         universe = bg.pid.tran$ENTREZID)
      
      # p1=dotplot(ego_BP)
      p1=barplot(ego_BP, showCategory=items, color="pvalue")
      
      return(list(barplot=p1,GO_obj=ego_BP,GO_table=summary(ego_BP)))
      
    }
    
    output_name = "./"
    
    ## main figure
    {
      # Cd4+Cd25+Klrg1- vs Cd4+Cd25+Klrg1+
      comparision_name = "Fig6_F_Cd4+Cd25+Klrg1-vsCd4+Cd25+Klrg1+"
      ## Cd25+Klrg1- vs Cd25+Klrg1+ down
      deg.df = read.csv(file="Fig6_E_Cd4+Cd25+Klrg1-vsCd4+Cd25+Klrg1+_2022-10-10.csv",header = T,row.names = 1)
      deg.down.df = deg.df %>% filter(P.Value<=0.05 & logFC<=-1)
      GO_BP_Outout.obj = run_GO_BP_BG(row.names(deg.down.df),row.names(deg.df))
      
      pdf(file=paste0(output_name,"/",comparision_name,"_Sig_Down_GO-BP_Plot_",Sys.Date(),".pdf"),
          width = 7,height = 7)  
      print(GO_BP_Outout.obj$barplot)
      print(plotGOgraph(GO_BP_Outout.obj$GO_obj))
      dev.off()
      
      write.csv(GO_BP_Outout.obj$GO_table,paste0(output_name,"/",comparision_name,"_Sig_Down_GO-BP_",Sys.Date(),".csv"),
                row.names =FALSE)
      
      ## Cd25+Klrg1- vs Cd25+Klrg1+ up
      deg.up.df = deg.df %>% filter(P.Value<=0.05 & logFC>=1)
      GO_BP_Outout.obj = run_GO_BP_BG(row.names(deg.up.df),row.names(deg.df))
      
      pdf(file=paste0(output_name,"/",comparision_name,"_Sig_Up_GO-BP_Plot_",Sys.Date(),".pdf"),
          width = 7,height = 7)  
      print(GO_BP_Outout.obj$barplot)
      print(plotGOgraph(GO_BP_Outout.obj$GO_obj))
      dev.off()
      
      write.csv(GO_BP_Outout.obj$GO_table,paste0(output_name,"/",comparision_name,"_Sig_Up_GO-BP_",Sys.Date(),".csv"),
                row.names =FALSE)
      
    }
    
    # supp figure
    {
      # Cd4+Cd25+ vs Cd4+Klrg1+
      {
        comparision_name = "Cd4+Cd25+vsCd4+Klrg1+"
        ## Cd25+Klrg1- vs Cd25+Klrg1+ down
        deg.df = read.csv(file="2_DEG_Cd25vsKlrg1_All/Cd4+Cd25+vsCd4+Klrg1+_2022-07-27.csv",header = T,row.names = 1)
        deg.down.df = deg.df %>% filter(P.Value<=0.05 & logFC<=-1)
        GO_BP_Outout.obj = run_GO_BP_BG(row.names(deg.down.df),row.names(deg.df))
        
        pdf(file=paste0(output_name,"/",comparision_name,"_Sig_Down_GO-BP_Plot_",Sys.Date(),".pdf"),
            width = 7,height = 7)  
        print(GO_BP_Outout.obj$barplot)
        print(plotGOgraph(GO_BP_Outout.obj$GO_obj))
        dev.off()
        
        write.csv(GO_BP_Outout.obj$GO_table,paste0(output_name,"/",comparision_name,"_Sig_Down_GO-BP_",Sys.Date(),".csv"),
                  row.names =FALSE)
        
        ## Cd25+Klrg1- vs Cd25+Klrg1+ up
        deg.up.df = deg.df %>% filter(P.Value<=0.05 & logFC>=1)
        GO_BP_Outout.obj = run_GO_BP_BG(row.names(deg.up.df),row.names(deg.df))
        
        pdf(file=paste0(output_name,"/",comparision_name,"_Sig_Up_GO-BP_Plot_",Sys.Date(),".pdf"),
            width = 7,height = 7)  
        print(GO_BP_Outout.obj$barplot)
        print(plotGOgraph(GO_BP_Outout.obj$GO_obj))
        dev.off()
        
        write.csv(GO_BP_Outout.obj$GO_table,paste0(output_name,"/",comparision_name,"_Sig_Up_GO-BP_",Sys.Date(),".csv"),
                  row.names =FALSE)
        
      }
      
      # Cd4+Klrg1- vs Cd4+Klrg1+
      {
        comparision_name = "Cd4+Klrg1-vsCd4+Klrg1+"
        #down
        deg.df = read.csv(file="2_DEG_Cd25vsKlrg1_All/Cd4+Klrg1-vsCd4+Klrg1+_2022-07-27.csv",header = T,row.names = 1)
        deg.down.df = deg.df %>% filter(P.Value<=0.05 & logFC<=-1)
        GO_BP_Outout.obj = run_GO_BP_BG(row.names(deg.down.df),row.names(deg.df))
        
        pdf(file=paste0(output_name,"/",comparision_name,"_Sig_Down_GO-BP_Plot_",Sys.Date(),".pdf"),
            width = 7,height = 7)  
        print(GO_BP_Outout.obj$barplot)
        print(plotGOgraph(GO_BP_Outout.obj$GO_obj))
        dev.off()
        
        write.csv(GO_BP_Outout.obj$GO_table,paste0(output_name,"/",comparision_name,"_Sig_Down_GO-BP_",Sys.Date(),".csv"),
                  row.names =FALSE)
        
        #up
        deg.up.df = deg.df %>% filter(P.Value<=0.05 & logFC>=1)
        GO_BP_Outout.obj = run_GO_BP_BG(row.names(deg.up.df),row.names(deg.df))
        
        pdf(file=paste0(output_name,"/",comparision_name,"_Sig_Up_GO-BP_Plot_",Sys.Date(),".pdf"),
            width = 7,height = 7)  
        print(GO_BP_Outout.obj$barplot)
        print(plotGOgraph(GO_BP_Outout.obj$GO_obj))
        dev.off()
        
        write.csv(GO_BP_Outout.obj$GO_table,paste0(output_name,"/",comparision_name,"_Sig_Up_GO-BP_",Sys.Date(),".csv"),
                  row.names =FALSE)
        
      }
      
      # Cd8+Klrg1- vs Cd8+Klrg1+
      {
        comparision_name = "Cd8+Klrg1-vsCd8+Klrg1+"
        #down
        deg.df = read.csv(file="2_DEG_Cd25vsKlrg1_All/Cd8+Klrg1-vsCd8+Klrg1+_2022-07-27.csv",header = T,row.names = 1)
        deg.down.df = deg.df %>% filter(P.Value<=0.05 & logFC<=-1)
        GO_BP_Outout.obj = run_GO_BP_BG(row.names(deg.down.df),row.names(deg.df))
        
        pdf(file=paste0(output_name,"/",comparision_name,"_Sig_Down_GO-BP_Plot_",Sys.Date(),".pdf"),
            width = 7,height = 7)  
        print(GO_BP_Outout.obj$barplot)
        print(plotGOgraph(GO_BP_Outout.obj$GO_obj))
        dev.off()
        
        write.csv(GO_BP_Outout.obj$GO_table,paste0(output_name,"/",comparision_name,"_Sig_Down_GO-BP_",Sys.Date(),".csv"),
                  row.names =FALSE)
        
        #up
        deg.up.df = deg.df %>% filter(P.Value<=0.05 & logFC>=1)
        GO_BP_Outout.obj = run_GO_BP_BG(row.names(deg.up.df),row.names(deg.df))
        
        pdf(file=paste0(output_name,"/",comparision_name,"_Sig_Up_GO-BP_Plot_",Sys.Date(),".pdf"),
            width = 7,height = 7)  
        print(GO_BP_Outout.obj$barplot)
        print(plotGOgraph(GO_BP_Outout.obj$GO_obj))
        dev.off()
        
        write.csv(GO_BP_Outout.obj$GO_table,paste0(output_name,"/",comparision_name,"_Sig_Up_GO-BP_",Sys.Date(),".csv"),
                  row.names =FALSE)
        
      }
    }
    
    # main, comparison go analysis
    {
      all.df = read.csv(file="Fig6_E_Cd4+Cd25+Klrg1-vsCd4+Cd25+Klrg1+_2022-10-10.csv",header = T,row.names = 1)
      deg.df = all.df
      deg.df$group = ""
      deg.df$group[deg.df$P.Value<=0.05 & deg.df$logFC<=-1] = "Cd4+Cd25+Klrg1-"
      deg.df$group[deg.df$P.Value<=0.05 & deg.df$logFC>=1] = "Cd4+Cd25+Klrg1+"
      deg.df = deg.df %>% filter(group !="")
      
      sig.pid.tran <- bitr(row.names(deg.df), fromType = "UNIPROT", toType = c("ENSEMBL", "ENTREZID","SYMBOL"), 
                           OrgDb = org.Mm.eg.db)
      bg.pid.tran <- bitr(row.names(all.df), fromType = "UNIPROT", toType = c("ENSEMBL", "ENTREZID","SYMBOL"), 
                          OrgDb = org.Mm.eg.db)
      
      deg.df = deg.df[sig.pid.tran$UNIPROT,]
      deg.df$Entrez = sig.pid.tran$ENTREZID
      formula_res <- compareCluster(Entrez~group, data=deg.df, fun="enrichGO", 
                                    OrgDb = org.Mm.eg.db,
                                    ont = "BP", pAdjustMethod = "BH",
                                    pvalueCutoff = 1,qvalueCutoff = 1,
                                    universe = bg.pid.tran$ENTREZID
      )
      
      head(formula_res)
      pdf(file="Fig6_F_Cd4+Cd25+Klrg1-vsCd4+Cd25+Klrg1+_Sig_Comparison_GO-BP_DotPlot_2022-07-28.pdf",
          width=11,height=8)
      dotplot(formula_res,color="qvalue",showCategory = 10,label_format=50,font.size=16)
      dev.off()
    }
    
    
    
  }
  

}
###-----------------------------------------------------

### estimate copy num
{
  Obj.list = readRDS(file="Validation_Cd25_Klrg1_All_Obj_2022-10-10.rds")
  ##> for copynum
  temp.df = Obj.list$identified
  mmu_mass.df = read.table(file="../../resource/uniprot-download_true_fields_accession_2Cmass_format_tsv_query__28mo-2022.10.10-07.07.18.89.tsv",header = T)
  temp.df$pid = row.names(temp.df)
  temp_anno.df = temp.df %>% merge(mmu_mass.df,by.x="pid",by.y="Entry",all.x=F,all.y=F)
  write.csv(x = temp_anno.df,
            file=paste("Identified_Intentisty_Mass_",dim(temp_anno.df)[1],"_",Sys.Date(),".csv",sep=""),quote = F,row.names = F)
  
   
}
  
###---------------------------
### Not use
{
  ### 04-Expression Distribution
  {
    {
      rm(list=ls())
      Obj.list = readRDS(file = "Validation_Cd25_Klrg1_All_Obj_2022-10-10.rds")
      
      temp.df = Obj.list$mann_imputate_log2 %>% 
        gather(key="name",value="value")
      
      temp.df = temp.df %>% merge(Obj.list$meta,by.x="name",by.y="sample_id")
      temp.df$name = factor(temp.df$name,levels=Obj.list$meta$sample_id)
      
      temp.p1=ggplot(data = temp.df,mapping = aes(x=name,y=value,fill=group_color)) + 
        geom_boxplot(fill=Obj.list$meta$group_color,outlier.colour = "black",color="black",
                     outlier.size = 0.1) +
        # geom_jitter(color="black", size=1, alpha=1) +
        theme_bw()+
        theme(
          legend.position="none",
          plot.title = element_text(size=18,hjust = 0.5),axis.text = element_text(size=10),
          axis.title = element_text(size=16),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=10),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          legend.key=element_blank()
        ) +
        stat_boxplot(geom = "errorbar", width = 0.15) +
        labs(title=paste("Expression distribution of all samples"),x="Cell types", y = "Log2(LFQ.Intensity")
      
      
      pdf(file = paste("1_ProteinExpression_Cd25vsKlrg1_All/5_Quantified_PG_Expression_Distribution/ProteinExpressionDistribution_MannImputated_Boxplot_",Sys.Date(),".pdf",sep=""),
          width=8,height=6)
      print(temp.p1)
      dev.off()
    }
    
  }
  
  ### 05-Correlation of protein
  {
    ## V1:no imputated
    {
      rm(list=ls())
      Obj.list = readRDS(file = "Validation_Cd25_Klrg1_All_Obj_2022-10-10.rds")
      
      temp.df = Obj.list$no_imputate_log2
      temp.cor = cor(temp.df[])
      temp.cor.anno = Obj.list$meta %>% select(group)
      
      # mycolors <- levels(Obj.list$meta$group_color)
      # names(mycolors) <- levels(temp.cor.anno$group)
      # mycolors <- list(group = mycolors)
      
      pheatmap(temp.cor,cluster_rows = F,cluster_cols = F,width = 7,height = 6,#cellwidth = 3,cellheight = 3,
               show_rownames = F,show_colnames = T,border_color = NA,
               annotation_col = temp.cor.anno,annotation_names_row = F,
               annotation_names_col = F,
               annotation_row = temp.cor.anno,#annotation_colors = mycolors,
               annotation_legend = F,
               # gaps_row = c(5, 24,44,69),gaps_col = c(5, 24,44,69),
               # legend_breaks = seq(0.6,1,by = 0.1),
               # legend_labels = c("0.6", "0.7", "0.8", "0.9", "Pearson Correlation\n"),
               # legend = TRUE, 
               # color = colorRampPalette(c("darkgrey","white","darkgoldenrod1","brown3"))(100),
               filename = paste("1_ProteinExpression_Cd25vsKlrg1_All/6_Quantified_PG_Expreesion_Corr/Corr_Pheatmap_NoImputate_",Sys.Date(),".pdf",sep="")
      )
    }
    
    ## V2:mann imputated
    {
      rm(list=ls())
      Obj.list = readRDS(file = "Validation_Cd25_Klrg1_All_Obj_2022-10-10.rds")
      
      temp.df = Obj.list$mann_imputate_log2
      temp.cor = cor(temp.df[])
      temp.cor.anno = Obj.list$meta %>% select(group)
      
      mycolors <- levels(Obj.list$meta$group_color)
      names(mycolors) <- levels(temp.cor.anno$group)
      mycolors <- list(group = mycolors)
      
      pheatmap(temp.cor,cluster_rows = F,cluster_cols = F,width = 7,height = 6,#cellwidth = 3,cellheight = 3,
               show_rownames = F,show_colnames = F,border_color = NA,
               annotation_col = temp.cor.anno,annotation_names_row = F,annotation_names_col = F,
               annotation_row = temp.cor.anno,annotation_colors = mycolors,
               annotation_legend = F,
               gaps_row = c(5, 24,44,69),gaps_col = c(5, 24,44,69),
               # legend_breaks = seq(0.6,1,by = 0.1),
               # legend_labels = c("0.6", "0.7", "0.8", "0.9", "Pearson Correlation\n"),
               # legend = TRUE, 
               color = colorRampPalette(c("darkgrey","white","darkgoldenrod1","brown3"))(100),
               filename = paste("6_Quantified_PG_Expreesion_Corr/Corr_Pheatmap_MannImputate_",Sys.Date(),".pdf",sep="")
      )
    }  
  }  
  
  ### 06-known marker 
  {
    ### Pheatmap,cell-specific marker
    {
      rm(list=ls())  
      Obj.list = readRDS(file = "Validation_Cd25_Klrg1_All_Obj_2022-10-10.rds")
      
      temp_PCC_marker = c("Epcam","Msln","Mmp7","Krt18")  
      temp_CAF_marker = c("Pdpn","Dcn",#pan-CAF
                          "Has1","Pdgfra","C8b","C3",#iCAF,
                          "Ly6c2","Tnc","Tagln","Ccn2","C2","Acta2","Thy1","Lrrc15",#myCAF, Ccn2,alias:Ctgf;aSMA,alias:Acta2
                          "Cd74"#,"H2-Aa","H2-Ab1" #apCAF,Cd74,alias MHCII;H2-Aa,P14434;H2-Ab1 P14483
      )
      temp_T_marker = c("Ptprc",#pan-immune,alias:Cd45
                        "Cd3e","Cd3g",#T-cell
                        "Cd4",
                        "Cd8a","Cd8b",
                        "Foxp3","Il2ra"#Il2ra, alias Cd25
      )
      temp_B_marker= c("Cd19","Cd79a","Cd79b","Ms4a1")
      temp_Myeloid_marker=c("Itgam",#Myeloid alias Cd11b
                            "Adgre1","Fn1","Apoe","Prg4","C1qb","C1qc","Saa3","Serpinb2",#Macrophage Gpf480 alias:F4/80,
                            "Ly6g","Cd14","Fcgr3","Fcgr4","Ccr2",#Mono
                            "Itgax","Slc46a3","Gm2a","Adam23","Naaa","Anpep","Pip4k2a",#DCs,Itgax alias:Cd11c
                            "S100a8","Ngp","S100a11","S100a9","Chi3l1","Lcn2","Mmp9","Pglyrp1"#Neu
      )
      
      temp.used.marker = c(temp_PCC_marker,
                           temp_CAF_marker,
                           temp_T_marker,
                           temp_B_marker,
                           temp_Myeloid_marker)
      temp.used.marker.df = Obj.list$annotation %>% filter(genename %in% temp.used.marker)
      row.names(temp.used.marker.df) = temp.used.marker.df$genename
      temp.used.marker.df = temp.used.marker.df[temp.used.marker,]
      
      temp.pheatmap.df = Obj.list$min_10_imputate_log2_row_means[temp.used.marker.df$pid,]
      row.names(temp.pheatmap.df) = temp.used.marker
      
      
      temp.samples.anno = Obj.list$meta %>% select(group) %>% unique()
      row.names(temp.samples.anno) = temp.samples.anno$group
      # temp.samples.anno$group = NULL
      
      mycolors <- levels(Obj.list$meta$group)
      names(mycolors) <- levels(temp.samples.anno$group)
      mycolors <- list(group = mycolors)
      
      pheatmap::pheatmap(temp.pheatmap.df,cellwidth = 40,cellheight = 15,
                         scale = "none",
                         cluster_rows=FALSE,cluster_cols = FALSE,angle_col=45,
                         annotation_names_row = F,annotation_names_col = F,
                         # gaps_row = c(4, 13,21,25),#border_color = NA,
                         # gaps_col = c(1, 5,10,15),
                         # width = 50,height = 400,
                         fontsize = 16,display_numbers = FALSE,
                         # annotation_col = temp.samples.anno,
                         # annotation_colors = mycolors,
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         filename = paste("1_ProteinExpression_Cd25vsKlrg1_All/7_Quantified_Known_Markers_Pheatmap/Known_Markers_All_Pheatmap_",Sys.Date(),".pdf",sep=""),
                         width = 15,height = 15
                         
      )
      
    }
    
    ### Pheatmap,select marker
    {
      rm(list=ls())  
      Obj.list = readRDS(file = "Validation_Cd25_Klrg1_All_Obj_2022-10-10.rds")
      # temp_PCC_marker = c("Epcam")  
      # temp_CAF_marker = c("Dcn",#pan-CAF
      #                     "Has1",#iCAF,
      #                     "Tagln",#myCAF, Ccn2,alias:Ctgf;aSMA,alias:Acta2
      #                     "H2-Aa" #apCAF,Cd74,alias MHCII;H2-Aa,P14434;H2-Ab1 P14483
      # )
      temp_T_marker = c("Ptprc",#pan-immune,alias:Cd45
                        "Cd3e",#"Cd3g",#T-cell
                        "Cd4",
                        "Cd8a","Cd8b",
                        "Foxp3","Klrg1","Tnfrsf18"
      )
      # temp_B_marker= c("Cd19")
      # temp_Myeloid_marker=c("Itgam",#Myeloid alias Cd11b
      #                       "S100a8",#Neu
      #                       "Cd14",#Mono.Fcgr3,alias Cd16;Fcgr4:Cd16a
      #                       "Adgre1",#Macrophage Gpf480 alias:F4/80,
      #                       "Itgax"#DCs,Itgax alias:Cd11c
      #                       
      # )
      
      temp.used.marker = c(
        temp_T_marker
      )
      temp.used.marker.df = Obj.list$annotation %>% filter(genename %in% temp.used.marker)
      # temp.used.marker.df = temp.used.marker.df[-which(temp.used.marker.df$pid=="P14437"),]
      row.names(temp.used.marker.df) = temp.used.marker.df$genename
      temp.used.marker.df = temp.used.marker.df[temp.used.marker,]
      
      temp.pheatmap.df = Obj.list$min_10_imputate_log2_row_means[temp.used.marker.df$pid,]
      row.names(temp.pheatmap.df) = temp.used.marker.df$genename
      
      
      # temp.samples.anno = Obj.list$meta %>% select(group,cell_classification) %>% unique()
      # row.names(temp.samples.anno) = temp.samples.anno$group
      # temp.samples.anno$group = NULL
      
      # mycolors <- levels(Obj.list$meta$cell_classification_color)
      # names(mycolors) <- levels(temp.samples.anno$cell_classification)
      # mycolors <- list(cell_classification = mycolors)
      
      pheatmap::pheatmap(temp.pheatmap.df,cellwidth = 25,cellheight = 15,
                         scale = "none",border_color = "white",
                         cluster_rows=FALSE,cluster_cols = FALSE,angle_col=45,
                         annotation_names_row = F,annotation_names_col = F,
                         # gaps_row = c(1, 5,8,9),#border_color = NA,
                         # gaps_col = c(1, 5,8,9),
                         # width = 50,height = 400,
                         fontsize = 12,display_numbers = FALSE,
                         # annotation_col = temp.samples.anno,
                         # annotation_colors = mycolors,
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         filename = paste0("1_ProteinExpression_Cd25vsKlrg1_All/7_Quantified_Known_Markers_Pheatmap/Selected_Known_Markers_Pheatmap_",Sys.Date(),".pdf"),
                         width = 6,height = 5)
      
      
    }
    
  }
  ### 10-DEG cluster, not run
  {
    ## Cd25+Klrg1- vs Cd25+Klrg1+ downregulate, GO BP
    {
      rm(list=ls())
      deg = read.csv(file="2_DEG_Cd25vsKlrg1_All/KPC_425KNvsKPC_425KP_VolcanoPlot_2022-07-22.csv",header = T,row.names = 1)
      deg.down = deg %>% filter(P.Value<=0.05 & logFC>=1)
    
      Obj.list = readRDS("Validation_Cd25_Klrg1_All_Obj_2022-10-10.rds")
  
      meta = Obj.list$meta %>% filter(group %in% c("KPC_425KN","KPC_425KP"))
      
      data.df = Obj.list$min_10_imputate_log2[row.names(deg.down),meta$sample_id]
      row.names(data.df) = deg.down$gene
      # row.anno = data.frame(row.names =row.names(deg.down),genename=deg.down$gene)
      
      
      pheatmap::pheatmap(data.df,cellwidth = 40,cellheight = 15,
                         scale = "row",
                         cluster_rows=TRUE,cluster_cols = FALSE,angle_col=45,
                         annotation_names_row = F,annotation_names_col = F,
                         # gaps_row = c(4, 13,21,25),#border_color = NA,
                         # gaps_col = c(1, 5,10,15),
                         # width = 50,height = 400,
                         fontsize = 16,display_numbers = FALSE,
                         # annotation_row =  row.anno,
                         # annotation_colors = mycolors,
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                         # filename = paste("1_ProteinExpression_Cd25vsKlrg1_All/7_Quantified_Known_Markers_Pheatmap/Known_Markers_All_Pheatmap_",Sys.Date(),".pdf",sep=""),
                         width = 15,height = 15
                         
      )
      
    
      # Obj.list$`Cd4+Cd25+Klrg1-_up` = ego_BP
      
    }
    
  }

}

#------------------
## color

# library(RColorBrewer)
# display.brewer.all(type = "seq")
# display.brewer.pal(9, "YlOrRd")
# display.brewer.pal(9, "Blues")
# brewer.pal(9, "Blues")
# temp.df = temp.df[,meta[order(meta$group),"id"]]

#---------------------
# library(enrichplot)
# library(ggnewscale)
# xx <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db")
# xx2 <- pairwise_termsim(xx)
