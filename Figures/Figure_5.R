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

# Figure_5d, heat map showing the Sig. proteins based on one-vs-the rest strategy
{
    rm(list=ls())
    Obj.list = readRDS(file = "ct_dataset_01.rds")
    data <- Obj.list$Mini_impute_log2
    deg <- Obj.list$dep_df
    meta = Obj.list$meta
    
    plot_data.df = data[unique(deg$pid),]
    plot_data.median=data.frame(pid=row.names(plot_data.df))
    for(i in levels(meta$CellType)){
      temp_SampleID = meta[meta$CellType  %in% i,"SampleID"]
      plot_data.median[,i]=apply(plot_data.df[,temp_SampleID],1,mean)
    }
    plot_data.median$genename = Obj.list$annotation[Obj.list$annotation$pid %in% plot_data.median$pid,"genename"]
        
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
    pdf(file=paste0("Figure_5d_Pheatmap_Alldeg_",Sys.Date(),".pdf"))
    print(p1)
    dev.off()
}

# Figure_5e, Enrichment Analysis based on CellLineage Sig. proteins
{
  rm(list=ls())
  Obj.list = readRDS(file = "ct_dataset_01.rds")
    
  data <- Obj.list$dep_df
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
  write.csv(formula_res_cutoff@compareClusterResult,file=paste0("Figure_5e_Overlap_CellLineage_Specific_GO-BP_DotPlot_",Sys.Date(),".csv"))
    
}
    
    
