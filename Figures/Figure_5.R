library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)

library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("")
# Figure_5b, bar plots of identified and quantified proteins, respectively
{
  rm(list=ls())
  Obj.list = readRDS(file = "ct_dataset_01.rds")
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
            file = paste0("Fig5_b_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".csv"),
            fileEncoding = "UTF-8",row.names = F,quote = TRUE)
  
  pdf(paste0("Fig5_b_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".pdf"),
        width=6,height = 4)
  print(temp.p1)
  dev.off()
}    

# Figure_12c, heat map of the Pearson Corr. Coef of all samples
{
  rm(list=ls())
  Obj.list = readRDS(file = "ct_dataset_01.rds")
  
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

  # begin plot
  p1 = pheatmap(temp.cor,cluster_rows = F,cluster_cols = F,width = 7,height = 6,#cellwidth = 3,cellheight = 3,
                show_rownames = T,show_colnames = F,border_color = NA,
                annotation_col = temp.cor.anno,annotation_names_row = F,annotation_names_col = F,
                annotation_row = temp.cor.anno,annotation_colors = mycolors,
                annotation_legend = F,fontsize_row = 6,fontsize_col = 6,
                color = colorRampPalette(c("darkgrey","white","darkgoldenrod1","brown3"))(100),silent = TRUE)
  pdf(file = paste0("Fig5_D_Corr_Pheatmap_Min10_", Sys.Date(),".pdf"), width=7,heigh=6)
  print(p1)
  dev.off()
  
  write.csv(temp.cor, file=paste0("Fig5_D_Corr_Pheatmap_Min10_", Sys.Date(),".csv"))
  
}

# Pheatmap for visualization of known markers
{
  rm(list=ls())  
  Obj.list = readRDS(file = "ct_dataset_01.rds")
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

# PCA showing the distribution of four lineages
{  
  rm(list=ls())
  Obj.list = readRDS(file = "ct_dataset_01.rds")
  temp.df = Obj.list$min_10_imputate_log2 %>% t() %>% as.data.frame()
  temp.df.anno = Obj.list$meta
  
  pca_result <- prcomp(temp.df,scale=T,center = T)
  
  df1 = as.data.frame(pca_result$x)
  summ1 <- summary(pca_result)
  xlab1 <- paste0("PC1 (",round(summ1$importance[2,1]*100,2),"%)")
  ylab1 <- paste0("PC2 (",round(summ1$importance[2,2]*100,2),"%)")

  
  p1=ggplot(data = df1,aes(x = PC1,y = PC2,color = temp.df.anno$CellType))+
    stat_ellipse(aes(fill = temp.df.anno$CellLineage),
                 type = "norm",geom = "polygon",alpha = 0.25,
                 color = NA
    ) + 
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
  
  
  pdf(file=paste0("Fig5_C_Quantified_PCA_CellLineage_MinImputated_",Sys.Date(),".pdf"),
      width = 7,height=5)
  print(p1)
  dev.off()
  
  write.csv(x = df1,file=paste0("Fig5_C_Quantified_PCA_CellLineage_MinImputated_",Sys.Date(),".csv"))
}

# Pheatmap showing the dep based on one-vs-the rest strategy
{
    rm(list=ls())
    Obj.list = readRDS(file = "ct_dataset_01.rds")
    data <- Obj.list$min_10_imputate_log2
    deg <- Obj.list$deg_Fig3
    meta = Obj.list$meta
    
    plot_data.df = data[unique(deg$pid),]
    plot_data.median=data.frame(pid=row.names(plot_data.df))
    for(i in levels(meta$CellType)){
      temp_SampleID = meta[meta$CellType  %in% i,"SampleID"]
      plot_data.median[,i]=apply(plot_data.df[,temp_SampleID],1,mean)
    }
    plot_data.median$genename = Obj.list$annotation[Obj.list$annotation$pid %in% plot_data.median$pid,"genename"]
    
    write.csv(plot_data.median,file = paste0("Fig3_E_Pheatmap_plot_data_",Sys.Date(),".csv"))
    
    ## add annotation
    annotation_col = data.frame(CellType = meta$CellType)
    row.names(annotation_col) = meta$SampleID
    anno_colors = list(CellType = levels(meta$CellType_Color))
    names(anno_colors$CellType) = levels(meta$CellType)
        
    p1 = pheatmap(plot_data.df,cluster_rows = F,show_rownames = F,cluster_cols = F,
             scale = "row",border_color = NA,legend = F,cellwidth = 5,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             annotation_col = annotation_col,show_colnames = F,legend_labels = F,annotation_legend = F,labels_col = NA,labels_row = NA,
             annotation_names_col = FALSE,annotation_colors = anno_colors
    )
    pdf(file=paste0("Fig3_E_Pheatmap_Alldeg_",Sys.Date(),".pdf"))
    print(p1)
    dev.off()
  }

# Enrichment Analysis based CellLineage Sig. proteins
{
  rm(list=ls())
  Obj.list = readRDS(file = "ct_dataset_01.rds")
    
  data <- Obj.list$deg_Fig3
  meta <- Obj.list$meta
  
  PCC_specific = data %>% filter(CellType %in% c("PCC")) %>% select(pid) %>% unlist()
  CAF_specific = data %>% filter(CellType %in% c("CAF","iCAF","myCAF","apCAF")) %>% select(pid) %>% unlist()
  Lymphoid_specific = data %>% filter(CellType %in% c("T4","T8","Treg","B")) %>% select(pid) %>% unlist()
  Myeloid_specific = data %>% filter(CellType %in% c("MYE","NEU","DC","MAC","MO")) %>% select(pid) %>% unlist()
  
  data.df = data.frame(pid=c(PCC_specific,CAF_specific,Lymphoid_specific,Myeloid_specific),
                       CellLineage=c(rep("Cancer Cell",length(PCC_specific)),rep("Fibroblast",length(CAF_specific)),rep("Lymphoid Cell",length(Lymphoid_specific)),rep("Myeloid Cell",length(Myeloid_specific))))
  data.df$CellLineage = factor(data.df$CellLineage, levels=levels(meta$CellLineage))
  target_pid.df = data.df
  
  target.pid.tran <- bitr(target_pid.df$pid, fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  target_pid.df = target_pid.df %>% merge(target.pid.tran,by.x="pid",by.y="UNIPROT",all.x=F,all.y=F) 
  
  bg.df = Obj.list$no_imputate_log2
  bg.df = bg.df[(apply(bg.df,1,sum))!=0,]
  bg.pid.tran <- bitr(row.names(bg.df), fromType = "UNIPROT", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  
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

  # begion plot
  p1=dotplot(formula_res_cutoff,x="CellLineage",color="pvalue",showCategory = 3,label_format=35,font.size=12)   

  # set title color
  p2 = p1
  p2$theme$axis.text.x$angle=45
  p2$theme$axis.text.x$hjust = 0.5
  p2$theme$axis.text.x$vjust = 0.5
  p2$theme$axis.text.x$colour= levels(Obj.list$meta$CellLineage_Color)
  p2$labels$x=NULL
  pdf(file=paste0("Fig3_D_Overlap_CellLineage_Specific_GO-BP_DotPlot_",Sys.Date(),".pdf"),
      width=6,height=6)
  print(p2)
  dev.off()
  write.csv(formula_res_cutoff@compareClusterResult,file=paste0("Fig3_D_Overlap_CellLineage_Specific_GO-BP_DotPlot_",Sys.Date(),".csv"))
    
}
    
    
