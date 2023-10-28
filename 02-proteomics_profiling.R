'---
title: "02_ProteinProfiling"
author: "Li-Yuan"
date: "2022-11-28"
output: github_document
---
'


## Protein Profiling
## Loading packages and setting global params
library(dplyr)
library(ggplot2) 
library(RColorBrewer)

setwd("E:\\XYF（许燕芬）\\FACS_SISPROT\\workspace\\20221128_final\\02-proteomics_profiling/")
options(stringsAsFactors=F)
options(encoding = "UTF-8")

#Fig2,C 
## Bar plot showing the protein numbers of quantified and identified proteins 
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  temp.df = Obj.list$identified %>% t() %>% as.data.frame()
  
  meta.df = Obj.list$meta %>% select(SampleID,CellType) %>% unique()
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
      axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1,color = c(levels(Obj.list$meta$CellType_Color),"Black")),
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
    scale_y_continuous(expand = c(0,0),limits=c(0,8000),breaks=c(0,2000,4000,5000,6000,7000,8000),label=c(0,2000,4000,5000,6000,7000,8000))
    
  
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

#sFig2,C
## Box plot showing the protein numbers of quantified and identified proteins
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221010.rds")
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
    ggplot( aes(x=CellType, y=protein_num, fill=protein_loca)) +
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
                                 color = levels(Obj.list$meta$CellType_Color)),
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


#Fig2,D
library(pheatmap)
## Corr
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  
  temp.df = Obj.list$min_10_imputate_log2
  for(i in 1:dim(temp.df)[1]){
    temp.df[i,which(temp.df[i,]==0)]=NA
  }
  temp.cor = cor(temp.df[])
  # row_max = c()#0.8996916
  for(i in 1:dim(temp.cor)[1]){
    temp_row = temp.cor[i,]
    # row_max = c(row_max,max(temp_row[which(temp_row!=1)]))
    temp.cor[i,which(temp_row==1)]=0.8996916
    
  }
  temp.cor.anno = Obj.list$meta %>% select(CellType)
  row.names(temp.cor.anno) = Obj.list$meta$SampleID
  mycolors <- c(
    as.character(levels(Obj.list$meta$CellType_Color))
  )
  names(mycolors) <- levels(temp.cor.anno$CellType)
  mycolors <- list(CellType = mycolors)
  # names(temp.pheatmap.df) = c(paste("T8_",1:5,sep=""),paste("T4",1:5,sep=""),paste("Treg",1:5,sep=""),
  #                             paste("MAC",1:5,sep=""),paste("DC_",1:5,sep=""))
  
  # write.csv(temp.cor,file="02_ProteinProfiling/Corr/Corr_Pheatmap_table.csv")
  p1 = pheatmap(temp.cor,cluster_rows = F,cluster_cols = F,width = 7,height = 6,#cellwidth = 3,cellheight = 3,
                show_rownames = T,show_colnames = F,border_color = NA,
                annotation_col = temp.cor.anno,annotation_names_row = F,annotation_names_col = F,
                annotation_row = temp.cor.anno,annotation_colors = mycolors,
                annotation_legend = F,fontsize_row = 6,fontsize_col = 6,
                color = colorRampPalette(c("darkgrey","white","darkgoldenrod1","brown3"))(100),silent = TRUE)
  pdf(file = paste0("Fig2_D_Corr_Pheatmap_Min10_",
                    Sys.Date(),".pdf"), width=7,heigh=6)
  print(p1)
  dev.off()
  
  write.csv(temp.cor, file=paste0("Fig2_D_Corr_Pheatmap_Min10_",
                                  Sys.Date(),".csv"))
  
}

library(grid)
#Fig3,A
## Pheatmap
{
  rm(list=ls())  
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  temp_PCC_marker = c("Epcam")  
  temp_CAF_marker = c("Dcn",#pan-CAF
                      "Has1",#iCAF,
                      "Tagln",#myCAF, Ccn2,alias:Ctgf;aSMA,alias:Acta2
                      "H2-Aa" #apCAF,Cd74,alias MHCII;H2-Aa,P14434;H2-Ab1 P14483
  )
  temp_T_marker = c("Ptprc",#pan-immune,alias:Cd45
                    "Cd3e",#"Cd3g",#T-cell
                    "Cd4",
                    "Cd8a",
                    "Foxp3"
  )
  temp_B_marker= c("Cd19")
  temp_Myeloid_marker=c("Itgam",#Myeloid alias Cd11b
                        "S100a8",#Neu
                        "Cd14",#Mono.Fcgr3,alias Cd16;Fcgr4:Cd16a
                        "Adgre1",#Macrophage Gpf480 alias:F4/80,
                        "Itgax"#DCs,Itgax alias:Cd11c
                        
  )
  
  known_markers = c(temp_PCC_marker,
                    temp_CAF_marker,
                    temp_T_marker,
                    temp_B_marker,
                    temp_Myeloid_marker)
  known_markers.df = Obj.list$annotation %>% filter(genename %in% known_markers)
  known_markers.df = known_markers.df[-which(known_markers.df$pid=="P14437"),]
  row.names(known_markers.df) = known_markers.df$genename
  known_markers.df = known_markers.df[known_markers,] ## sort
  known_markers.df$commonly_genename = c("Epcam","Dcn","Has1","Tagln","H2-Aa","CD45","CD3","CD4","CD8a","Foxp3","CD19","CD11b","S100a8","CD14","F4/80","CD11c")
  
  known_markers_lfq.df = Obj.list$min_10_imputate_log2[known_markers.df$pid,]## all samples,markers
  
  ## row means
  meta.df = Obj.list$meta 
  known_marker_lfq_mean.df = data.frame(row.names = row.names(known_markers_lfq.df))
  for(i in levels(meta.df$CellType)){
    temp_index = which(meta.df$CellType==i)  
    known_marker_lfq_mean.df[,i] = rowMeans(known_markers_lfq.df[,temp_index],
                                            na.rm = TRUE)
  }
  
  row.names(known_marker_lfq_mean.df) = known_markers.df$commonly_genename
  
  p1=pheatmap::pheatmap(known_marker_lfq_mean.df,cellwidth = 18,cellheight = 15,
                        scale = "row",border_color = "white",
                        cluster_rows=FALSE,cluster_cols = FALSE,angle_col=45,
                        annotation_names_row = F,annotation_names_col = F,legend = T,
                        fontsize = 12,display_numbers = FALSE,
                        silent = TRUE,
                        color = colorRampPalette(c("navy", "white", "firebrick3"))(100), width = 9,height = 5)
  
  my_gtable = p1$gtable
  my_gtable$grobs[[2]]$gp=gpar(col=levels(Obj.list$meta$CellType_Color),fontsize=12)
  p1$gtable = my_gtable
  
  pdf(file = paste0("Fig3_A_Selected_Known_Markers_Pheatmap_",Sys.Date(),".pdf"),
      width=7,heigh=6)
  print(p1)
  dev.off()
  
  write.csv(known_marker_lfq_mean.df, file=paste0("Fig3_A_Selected_Known_Markers_Pheatmap_",Sys.Date(),".csv",sep=""),row.names = T)
  
}

library(ggrepel)
#Fig3,B
## Dynamic Range with known markers
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  
  temp.df = Obj.list$quantified
  temp.meta = Obj.list$meta
  temp.plot.data = matrix(nrow = length(levels(temp.meta$CellType)),ncol=length(temp.df[,1]),data = 0) %>% as.data.frame()
  names(temp.plot.data) = row.names(temp.df)
  row.names(temp.plot.data) = levels(temp.meta$CellType)
  
  temp.df = temp.df[,temp.meta$SampleID]
  for(i in levels(temp.meta$CellType)){
    temp_index = which(temp.meta$CellType==i)  
    temp_df = temp.df[,temp_index]
    temp.plot.data[i,] = rowMeans(temp_df,na.rm = TRUE)
    
  }
  temp.plot.data.log10 = temp.plot.data %>% t() %>% as.data.frame()
  
  temp.line.data = data.frame()
  for(i in (1:dim(temp.plot.data.log10)[2])){
    temp_df = data.frame(name=names(temp.plot.data.log10)[i],value=log10(temp.plot.data.log10[,i]+1),
                         pid=row.names(temp.plot.data.log10))
    temp_df = temp_df %>% filter(value>0) %>% mutate(rank=min_rank(desc(value)))
    temp.line.data = rbind(temp.line.data,temp_df)
    
  }
  temp.line.data$name = factor(temp.line.data$name, levels=levels(Obj.list$meta$CellType))
  
  temp.line.anno.data = temp.line.data %>% filter(
    (name=="PCC"&pid=="Q99JW5")|
      (name=="CAF"&pid=="P97321")|
      (name=="iCAF"&pid=="Q61647")|
      (name=="myCAF"&pid=="Q80YX1")|
      (name=="apCAF"&pid=="P04441")|
      
      (name=="B"&pid=="P25918")|
      (name=="T4"&pid=="P06332")|
      (name=="T8"&pid=="P01731")|
      (name=="Treg"&pid=="Q99JB6")|#Foxp3
      # (name=="NK" &pid=="Q08481")|
      
      (name=="MYE"&pid=="P05555")|#Cd11b,Itgam
      (name=="MAC"&pid=="Q61549")|#Adgre1,F4/80
      (name=="NEU"&pid=="P10810")|
      (name=="DC"&pid=="Q9QXH4")#Itgax,Cd11c
  )
  temp.line.anno.data$marker = c("Epcam","Fap","Has1","Tnc","Cd74",
                                 "Cd19","Cd4","Cd8a","Foxp3",#"Cd31",
                                 "Cd11b","F4/80","Cd14","Cd11c")
  set.seed(12311)
  temp.line.plot.with.anno = temp.line.data %>%
    ggplot( aes(x=rank, y=value, color=name, group=name)) +
    geom_line(size=1.5) +
    geom_text_repel( data=temp.line.anno.data, aes(x=rank,y=value,color=name, label = marker), hjust = -3, vjust = -1, size = 5.0,show.legend = F,min.segment.length = 1 ) +
    geom_point(data=temp.line.anno.data,aes(x=rank,y=value,fill=name,size=5.0),show.legend = F)+
    theme_classic()+
    theme(
      axis.title = element_text(size = 16),axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
      legend.text = element_text(size = 12),
      legend.key.size = unit(5,"mm"),
      # legend.spacing.y=unit(15,"lines"),
      legend.position = c(.94,.55),
      legend.background = element_blank(),panel.grid = element_blank(),
      plot.background = element_blank(),strip.background = element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank()
      
    )+
    guides(colour=guide_legend(title=""))+ #legend title
    scale_color_manual(values = levels(temp.meta$CellType_Color))+#group color
    scale_y_continuous(limits=c(2,8),breaks=c(2,4,6,8),label=c(2,4,6,8))+
    scale_x_continuous(limits=c(0,6000),breaks=c(0,2000,4000,6000,8000),label=c(0,2000,4000,6000,8000))+
    
    
    labs(x="Protein Abundance(rank)",y="Log10(LFQ.Intensity)")
  
  
  pdf(file=paste("Fig3_B_PG_Abundance_LinePlot_",Sys.Date(),".pdf",sep="")
      ,width = 6,height = 4)
  print(temp.line.plot.with.anno);
  dev.off()
  
  write.csv(temp.line.data, file=paste0("Fig3_B_PG_Abundance_LinePlot_",Sys.Date(),".csv",sep=""),row.names = F)
  write.csv(temp.line.anno.data, file=paste0("Fig3_B_PG_Abundance_LinePlot_Label_",Sys.Date(),".csv",sep=""),row.names = F)
  # print(temp.line.plot.with.anno)
}

#Fig3,C
## PCA showing the distribution of four lineages
{  
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  temp.df = Obj.list$min_10_imputate_log2 %>% t() %>% as.data.frame()
  temp.df.anno = Obj.list$meta
  
  pca_result <- prcomp(temp.df,scale=T,center = T)
  
  df1 = as.data.frame(pca_result$x)
  summ1 <- summary(pca_result)
  xlab1 <- paste0("PC1 (",round(summ1$importance[2,1]*100,2),"%)")
  ylab1 <- paste0("PC2 (",round(summ1$importance[2,2]*100,2),"%)")
  #
  p1=ggplot(data = df1,aes(x = PC1,y = PC2,color = temp.df.anno$CellType))+
    stat_ellipse(aes(fill = temp.df.anno$CellLineage),
                 type = "norm",geom = "polygon",alpha = 0.25,
                 color = NA
    )+ 
    geom_point(size = 3.5)+
    labs(x = xlab1,y = ylab1,color = "",title = "")+
    guides(fill = "none")+
    theme_bw()+
    scale_fill_manual(values = levels(temp.df.anno$CellLineage_Color))+
    scale_colour_manual(values = levels(temp.df.anno$CellType_Color))+
    theme(plot.background = element_blank(),legend.background = element_blank(),
          panel.background = element_blank(),panel.grid = element_blank(),
          axis.text = element_text(size = 12),axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
  
  
  pdf(file=paste0("Fig3_C_Quantified_PCA_CellLineage_MinImputated_",Sys.Date(),".pdf"),
      width = 7,height=5)
  print(p1)
  dev.off()
  
  write.csv(x = df1,file=paste0("Fig3_C_Quantified_PCA_CellLineage_MinImputated_",Sys.Date(),".csv"))
}

library(limma)
#Fig3,E
## Pheatmap with deg
{
  {
    
    rm(list=ls())
    FACS.MSFragger.pg.df = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
    
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
        
        # print(temp.df.select[c("P25918"),])
        
        temp.df.select.sig = temp.df.select %>% filter(P.Value<0.05 & logFC < -1)
        # temp.df.select.sig = temp.df.select.sig[order(abs(temp.df.select.sig$logFC),decreasing = TRUE),]
        
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
        # print(length(x))
      }
      temp.df.sig$logFC = abs(temp.logFC)
      temp.df.sig = temp.df.sig[order(temp.df.sig$logFC,decreasing = T),]
      temp.df.sig$CellType=levels(temp_meta$CellType)[1]
      return(temp.df.sig)
    }
    
    ### PCC vs Rest
    meta <- FACS.MSFragger.pg.df$meta
    
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
      
      temp_subdf = FACS.MSFragger.pg.df$min_10_imputate_log2 %>% 
        select(temp_submeta$SampleID)
      temp.out = run_OnevsRest(temp_submeta, temp_subdf)
      if(temp.out!=0){
        temp_out.df = rbind(temp_out.df,temp.out)  
      }
      
      print(i)
      # break
    }
    
    temp_out.df$genename = FACS.MSFragger.pg.df$annotation[temp_out.df$pid,"genename"]
    temp_out.df[temp_out.df$pid %in% FACS.MSFragger.pg.df$uniprot_mouse_loci$Entry,"Loca"]="Surface"
    write.csv(temp_out.df,file=paste0("OnevsRest_DEG_Min10Imputate",dim(temp_out.df)[1],"_",Sys.Date(),".csv"),row.names = F,na = "")
    
    # temp_out.df$CellType = factor(temp_out.df$CellType,levels=levels(FACS.MSFragger.pg.df$meta$CellType))
    # temp_out.df = temp_out.df[order(temp_out.df$CellType),]
    
    FACS.MSFragger.pg.df$deg_Fig3 = temp_out.df
    saveRDS(FACS.MSFragger.pg.df,file = "../FACS.MSFragger.Obj.20221128.rds")
  }
  
  {
    # Pheatmap
    rm(list=ls())
    Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221010.rds")
    data <- Obj.list$mann_imputate_log2
    deg <- Obj.list$deg_Fig3
    meta = Obj.list$meta
    
    plot_data.df = data[unique(deg$pid),]
    # meta <- Obj.list$meta %>% filter(CellType %in% c("PCC","CAF","iCAF","myCAF","apCAF"))  
    # deg <- Obj.list$deg%>% filter(CellType %in% c("PCC","CAF","iCAF","myCAF","apCAF")) 
    # # deg.pid = names(table(deg$pid))[which(table(deg$pid) ==1)]
    # 
    # plot.data = data[unique(deg.pid),meta$SampleID]
    # 
    plot_data.median=data.frame(pid=row.names(plot_data.df))
    for(i in levels(meta$CellType)){
      temp_SampleID = meta[meta$CellType  %in% i,"SampleID"]
      plot_data.median[,i]=apply(plot_data.df[,temp_SampleID],1,mean)
    }
    plot_data.median$genename = Obj.list$annotation[Obj.list$annotation$pid %in% plot_data.median$pid,"genename"]
    
    write.csv(plot_data.median,file = paste0("Fig3_E_Pheatmap_plot_data_",Sys.Date(),".csv"))
    
    ## annotation
    annotation_col = data.frame(
      CellType = meta$CellType
    )
    row.names(annotation_col) = meta$SampleID
    anno_colors = list(
      CellType = levels(meta$CellType_Color)
    )
    names(anno_colors$CellType) = levels(meta$CellType)
    
    pdf(file=paste0("Fig3_E_Pheatmap_Alldeg_",Sys.Date(),".pdf"))
    pheatmap(plot_data.df,cluster_rows = F,show_rownames = F,cluster_cols = F,
             scale = "row",border_color = NA,legend = F,cellwidth = 5,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(1000),
             annotation_col = annotation_col,show_colnames = F,legend_labels = F,annotation_legend = F,labels_col = NA,labels_row = NA,
             annotation_names_col = FALSE,annotation_colors = anno_colors
    )
    dev.off()
  }
}

#Fig3,D
## CellLineage Protein annotation
{
  library(DOSE)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221010.rds")
  data <- Obj.list$quantified
  meta <- Obj.list$meta
  
  PCC_data = data %>% select(meta %>% filter(CellLineage=="Cancer Cell") %>% select(SampleID) %>% unlist()) %>% rowSums()
  CAF_data = data %>% select(meta %>% filter(CellLineage=="Fibroblast") %>% select(SampleID) %>% unlist())%>% rowSums()
  Lymphoid_data = data %>% select(meta %>% filter(CellLineage=="Lymphoid Cell") %>% select(SampleID) %>% unlist())%>% rowSums()
  Myeloid_data = data %>% select(meta %>% filter(CellLineage=="Myeloid Cell") %>% select(SampleID) %>% unlist())%>% rowSums()
  
  PCC_data = names(PCC_data[which(PCC_data!=0)])
  CAF_data = names(CAF_data[which(CAF_data!=0)])
  Lymphoid_data = names(Lymphoid_data[which(Lymphoid_data!=0)])
  Myeloid_data = names(Myeloid_data[which(Myeloid_data!=0)])
  
  PCC_specific = setdiff(PCC_data,unique(c(CAF_data,Lymphoid_data,Myeloid_data)))#115
  CAF_specific = setdiff(CAF_data,unique(c(PCC_data,Lymphoid_data,Myeloid_data)))#371
  
  Lymphoid_specific = setdiff(Lymphoid_data,unique(c(CAF_data,PCC_data,Myeloid_data)))#104
  Myeloid_specific = setdiff(Myeloid_data,unique(c(CAF_data,Lymphoid_data,PCC_data)))#332
  
  data.df = data.frame(pid=c(PCC_specific,CAF_specific,Lymphoid_specific,Myeloid_specific),CellLineage=c(rep("Cancer Cell",length(PCC_specific)),rep("Fibroblast",length(CAF_specific)),rep("Lymphoid Cell",length(Lymphoid_specific)),rep("Myeloid Cell",length(Myeloid_specific))))
  data.df$CellLineage = factor(data.df$CellLineage, levels=levels(meta$CellLineage))
  
  data.df$genename = Obj.list$annotation[data.df$pid,"genename"]
  pid1.tran <- bitr(data.df$genename, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  
  pid1.tran.df = data.df %>% merge(pid1.tran,by.x="genename",by.y="SYMBOL",all.x=FALSE,all.y=TRUE) ##883
  
  
  pid2.tran = bitr(data.df$pid, fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  pid2.tran.df = data.df %>% merge(pid2.tran,by.x="pid",by.y="UNIPROT",all.x=FALSE,all.y=TRUE)##868
  pid.tran.df = unique(rbind(pid1.tran.df,pid2.tran.df))#919
  # setdiff(data.df$pid,pid.tran.df$pid)
  
  bg.pid.tran <- bitr(row.names(data), fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  formula_res <- compareCluster(ENTREZID~CellLineage, data=pid.tran.df, fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,
                                universe = bg.pid.tran$ENTREZID,readable=TRUE
  )
  
  head(formula_res)
  formula_res_cutoff = formula_res
  formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<=0.05,]
  
  
  p1=dotplot(formula_res_cutoff,x="CellLineage",color="pvalue",showCategory = 3,label_format=35,font.size=12) 
  # facet_grid(~cell_classification)
  
  
  p2 = p1
  p2$theme$axis.text.x$angle=45
  p2$theme$axis.text.x$hjust = 0.5
  p2$theme$axis.text.x$vjust = 0.5
  p2$theme$axis.text.x$colour= levels(Obj.list$meta$CellLineage_Color)
  p2$labels$x=NULL
  pdf(file=paste0("Fig3_D_Overlap_CellLineage_Specific_GO-BP_DotPlot_",Sys.Date(),".pdf"),
      width=6,height=6)
  p2
  dev.off()
  write.csv(formula_res_cutoff@compareClusterResult,file=paste0("Fig3_D_Overlap_CellLineage_Specific_GO-BP_DotPlot_",Sys.Date(),".csv"))
  
  
}

## Fig3,D,DEG
{
  library(DOSE)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  data <- Obj.list$deg_Fig3
  meta <- Obj.list$meta
  
  # PCC_data = data %>% select(meta %>% filter(CellLineage=="Cancer Cell") %>% select(SampleID) %>% unlist()) %>% rowSums()
  # CAF_data = data %>% select(meta %>% filter(CellLineage=="Fibroblast") %>% select(SampleID) %>% unlist())%>% rowSums()
  # Lymphoid_data = data %>% select(meta %>% filter(CellLineage=="Lymphoid Cell") %>% select(SampleID) %>% unlist())%>% rowSums()
  # Myeloid_data = data %>% select(meta %>% filter(CellLineage=="Myeloid Cell") %>% select(SampleID) %>% unlist())%>% rowSums()
  # 
  # PCC_data = names(PCC_data[which(PCC_data!=0)])
  # CAF_data = names(CAF_data[which(CAF_data!=0)])
  # Lymphoid_data = names(Lymphoid_data[which(Lymphoid_data!=0)])
  # Myeloid_data = names(Myeloid_data[which(Myeloid_data!=0)])
  # 
  PCC_specific = data %>% filter(CellType %in% c("PCC")) %>% select(pid) %>% unlist()
  CAF_specific = data %>% filter(CellType %in% c("CAF","iCAF","myCAF","apCAF")) %>% select(pid) %>% unlist()
  Lymphoid_specific = data %>% filter(CellType %in% c("T4","T8","Treg","B")) %>% select(pid) %>% unlist()
  Myeloid_specific = data %>% filter(CellType %in% c("MYE","NEU","DC","MAC","MO")) %>% select(pid) %>% unlist()
  
  data.df = data.frame(pid=c(PCC_specific,CAF_specific,Lymphoid_specific,Myeloid_specific),
                       CellLineage=c(rep("Cancer Cell",length(PCC_specific)),rep("Fibroblast",length(CAF_specific)),rep("Lymphoid Cell",length(Lymphoid_specific)),rep("Myeloid Cell",length(Myeloid_specific))))
  data.df$CellLineage = factor(data.df$CellLineage, levels=levels(meta$CellLineage))
  target_pid.df = data.df
  
  target.pid.tran <- bitr(target_pid.df$pid, fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  target_pid.df =target_pid.df %>% merge(target.pid.tran,by.x="pid",by.y="UNIPROT",all.x=F,all.y=F) 
  
  bg.df = Obj.list$no_imputate_log2
  bg.df = bg.df[(apply(bg.df,1,sum))!=0,]
  # row.names(bg.df) = gsub(row.names(bg.df),pattern = ";.*$",replacement = "")
  bg.pid.tran <- bitr(row.names(bg.df), fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  
  write.csv(file=paste0("Fig3D_deg",Sys.Date(),".csv"),x=target_pid.df)
  write.csv(file=paste0("Fig3D_bg",Sys.Date(),".csv"),x=bg.df)
  
  formula_res <- compareCluster(ENTREZID~CellLineage, data=target_pid.df, fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,
                                universe = bg.pid.tran$ENTREZID,readable=TRUE
  )
  
  formula_res_cutoff = formula_res
  formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<=0.05,]
  
  ID1 = formula_res_cutoff@compareClusterResult %>% 
    filter (CellLineage == "Cancer Cell") %>% select(ID) %>% head(n=3)
  ID2 = formula_res_cutoff@compareClusterResult %>% 
    filter (CellLineage == "Fibroblast") %>% select(ID) %>% head(n=3)
  ID3 = formula_res_cutoff@compareClusterResult %>% 
    filter (CellLineage == "Lymphoid Cell") %>% select(ID) %>% head(n=3)
  ID4 = formula_res_cutoff@compareClusterResult %>% 
    filter (CellLineage == "Myeloid Cell") %>% select(ID) %>% head(n=2)
  formula_res_cutoff2 = formula_res_cutoff@compareClusterResult %>% filter(ID %in% c(ID1$ID,ID2$ID,ID3$ID,ID4$ID))
  formula_res_cutoff@compareClusterResult = formula_res_cutoff2
  p1=dotplot(formula_res_cutoff,x="CellLineage",color="pvalue",showCategory = 3,label_format=35,font.size=12) 
  
  p2 = p1
  p2$theme$axis.text.x$angle=45
  p2$theme$axis.text.x$hjust = 0.5
  p2$theme$axis.text.x$vjust = 0.5
  p2$theme$axis.text.x$colour= levels(Obj.list$meta$CellLineage_Color)
  p2$labels$x=NULL
  pdf(file=paste0("Fig3_D_Overlap_CellLineage_Specific_GO-BP_DotPlot_",Sys.Date(),".pdf"),
      width=6,height=6)
  p2
  dev.off()
  write.csv(formula_res_cutoff@compareClusterResult,file=paste0("Fig3_D_Overlap_CellLineage_Specific_GO-BP_DotPlot_",Sys.Date(),".csv"))
  
  ##> GOMF
  formula_res <- compareCluster(ENTREZID~CellLineage, data=target_pid.df, fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,
                                ont = "MF", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,
                                universe = bg.pid.tran$ENTREZID,readable=TRUE
  )
  
  formula_res_cutoff = formula_res
  formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<=0.05,]
  
  
  p1=dotplot(formula_res_cutoff,x="CellLineage",color="pvalue",showCategory = 3,label_format=35,font.size=12) 
  
  p2 = p1
  p2$theme$axis.text.x$angle=45
  p2$theme$axis.text.x$hjust = 0.5
  p2$theme$axis.text.x$vjust = 0.5
  p2$theme$axis.text.x$colour= levels(Obj.list$meta$CellLineage_Color)
  p2$labels$x=NULL
  pdf(file=paste0("Fig3_D_Overlap_CellLineage_Specific_GO-MF_DotPlot_",Sys.Date(),".pdf"),
      width=6,height=6)
  p2
  dev.off()
  write.csv(formula_res_cutoff@compareClusterResult,file=paste0("Fig3_D_Overlap_CellLineage_Specific_GO-MF_DotPlot_",Sys.Date(),".csv"))
  
}

library(VennDiagram)
#Figs3,D
## Venn
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221010.rds")
  data <- Obj.list$quantified
  meta <- Obj.list$meta
  
  PCC_data = data %>% select(meta %>% filter(CellLineage=="Cancer Cell") %>% select(SampleID) %>% unlist()) %>% rowSums()
  CAF_data = data %>% select(meta %>% filter(CellLineage=="Fibroblast") %>% select(SampleID) %>% unlist())%>% rowSums()
  Lymphoid_data = data %>% select(meta %>% filter(CellLineage=="Lymphoid Cell") %>% select(SampleID) %>% unlist())%>% rowSums()
  Myeloid_data = data %>% select(meta %>% filter(CellLineage=="Myeloid Cell") %>% select(SampleID) %>% unlist())%>% rowSums()
  
  
  PCC_data = PCC_data[which(PCC_data!=0)]
  CAF_data = CAF_data[which(CAF_data!=0)]
  Lymphoid_data = Lymphoid_data[which(Lymphoid_data!=0)]
  Myeloid_data = Myeloid_data[which(Myeloid_data!=0)]
  
  #Make the plot
  venn.plot=venn.diagram(
    x = list(
      "PCC"=names(PCC_data),
      "CAF"=names(CAF_data),
      "Lymphoid Cell"=names(Lymphoid_data),
      "Myeloid Cell"=names(Myeloid_data)
    ),
    # category.names = c("PCC ()" , "CAF (663)" , "Lymphoid Cell (471)","A"),
    filename = NULL,
    output = TRUE ,
    imagetype="tiff" ,
    height = 480 ,
    width = 480 ,
    resolution = 600,
    compression = "lzw",
    lwd = 1,
    col=levels(meta$CellLineage_Color),
    fill = c(alpha("#F768A1",0.3), alpha('#88419D',0.3), alpha('#E31A1C',0.3),alpha('#006D2C',0.3)),
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(0, 50, 100,200),
    cat.dist = c(0.055, 0.055, 0.085,0.085),
    cat.fontfamily = "sans",
    cat.col = levels(meta$CellLineage_Color)#,
    # rotation = 1
  )
  
  pdf(file=paste("FigS3_D_PG_Quantified_VennPlot_CellLineage_",Sys.Date(),".pdf",sep=""),width = 4,height = 5)  
  grid.draw(venn.plot)
  dev.off()
  
}

## Expression Distribution
## Not used
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221010.rds")

  temp.df = Obj.list$mann_imputate_log2 %>%
    gather(key="name",value="value")

  temp.df = temp.df %>% merge(Obj.list$meta,by.x="name",by.y="SampleID")
  temp.df$name = factor(temp.df$name,levels=Obj.list$meta$SampleID)

  temp.p1=ggplot(data = temp.df,mapping = aes(x=name,y=value,fill=CellType_Color)) + geom_boxplot(fill=Obj.list$meta$CellType_Color,outlier.colour = "black",color="black",outlier.size = 0.1) +
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

