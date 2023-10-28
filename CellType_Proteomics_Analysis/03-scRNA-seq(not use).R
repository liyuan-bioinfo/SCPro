#'@author LiYuan
#'@desc comparision with scRNA-seq
#'@version final-version


library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyr)

options(stringsAsFactors = F)
options(encoding = "UTF-8")

setwd("E:\\XYF（许燕芬）\\FACS_SISPROT\\workspace\\20221128_final\\03-scRNA-seq/")

#panel a, tsne
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2022-09-14.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  
  ##> scRNA-seq
  scRNA_meta.df = TotalTissue.harmony@meta.data
  # scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  # scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  # scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  # scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
  # scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
  
  # scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
  # scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(Obj.list$meta$CellType))
  scRNA_meta.df$Main_subset_labels2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[which(scRNA_meta.df$RNA_snn_res.0.5==7),"Main_subset_labels2"]="EMT_like"
  
  ##> plot
  scRNA_meta.df$Main_subset_labels2=factor(scRNA_meta.df$Main_subset_labels2,levels=(
    c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC",
      "Endo","Acinar","Pericyte","EMT_like","Other")))
  # levels(Obj.list$meta$CellType_Color)
  # levels(Obj.list$meta$CellType)
  # ))
  TotalTissue.harmony@meta.data = scRNA_meta.df
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels2
  
  
  pdf(file=paste0("C_scRNA-seq_UMAP_",Sys.Date(),".pdf"),width=7,height = 5)
  DimPlot(TotalTissue.harmony,label = FALSE,cols = c(
    "#F768A1","#9EBCDA","#8C96C6", "#8C6BB1", 
    "#FED976","#FEB24C", "#FD8D3C" ,"#FC4E2A", 
    "#CCECE6", "#99D8C9", "#66C2A4","#41AE76", "#238B45",
    brewer.pal(9, "Greys")[4:8]
  ))
  dev.off()
  
  # write.csv(TotalTissue.harmony@meta.data,file="reanalysis.csv")
  saveRDS(TotalTissue.harmony,file=paste0("GSE129455_All_Viable_with_subset_",Sys.Date(),".rda"))
}

#panel a, David Tuevson umap
{
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2022-10-09.rda")
  Idents(TotalTissue.harmony) = TotalTissue.harmony$RNA_snn_res.0.5
  
  # subcell.Seurat <- RunUMAP(object = TotalTissue.harmony, dims = 1:40, reduction = 'harmony')
  
  DimPlot(TotalTissue.harmony, reduction = "umap",label = T, pt.size=0.5)
  
  MAC=c("Apoe","Saa3","C1qc") #0,1,2,4
  MYE=c("Ear2","Lyz1") #6
  B = c("Cd79a","Ly6d","Ms4a1") #15
  NEU = c("S100a8","S100a9","G0s2") #3
  Fib = c("Dcn","Col1a1","Col3a1") #11
  T_NK = c("Ccl5","Cd3g","Gzma","Nkg7") #12
  Ductal = c("Tff1","Krt18","Krt8","Krt19","Krt7","Epcam")#5
  Endo = c("Igfbp7","Plvap","Cd34")#16
  DCs = c("Ccl5","Ccr7","Ccl22")#14
  Acinar = c("Ctrb1","Prss2","Try5")#17
  Peri = c("Igfbp7","Rgs5","Acta2")#18
  EMT_like = c("Cdkn2a","S100a6","Igfbp4","Sparc","Vim","Spp1")#7
  
  DotPlot(TotalTissue.harmony, features = MYE)
  VlnPlot(TotalTissue.harmony, features = DCs)
  current.cluster.ids = c(0:18)
  current.cluster.ids[c(1,2,3,5,10,11)]="MAC"
  current.cluster.ids[c(7)]="MYE"
  current.cluster.ids[c(16)]="B"
  current.cluster.ids[c(4)]="NEU"
  current.cluster.ids[c(12)]="CAF"
  current.cluster.ids[c(13)]="T_NK"
  current.cluster.ids[c(6)]="Ductal"
  current.cluster.ids[c(17)]="Endo"
  current.cluster.ids[c(9,14,15)]="DCs"
  current.cluster.ids[c(18)]="Acinar"
  current.cluster.ids[c(19)]="Peri"
  current.cluster.ids[c(8)]="EMT_like"
  # Idents(TotalTissue.harmony) = TotalTissue.harmony$RNA_snn_res.0.5
  
  names(x = current.cluster.ids) <- levels(x = TotalTissue.harmony)
  TotalTissue.harmony <- RenameIdents(object = TotalTissue.harmony, current.cluster.ids)
  TotalTissue.harmony[["David_Tuevson_Clusters"]] = Idents(object = TotalTissue.harmony)
  

  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$David_Tuevson_Clusters=factor(scRNA_meta.df$David_Tuevson_Clusters,levels=(
    c("Ductal","CAF","T_NK","B","MYE",
      "NEU","MAC","DCs","Endo","Acinar","Peri","EMT_like"
      )))
  TotalTissue.harmony@meta.data = scRNA_meta.df
  
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels2
  p1=DimPlot(TotalTissue.harmony,label = TRUE,raster = T,cols = c(
    "#F768A1","#9EBCDA","#8C96C6", "#8C6BB1", 
    "#FED976","#FEB24C", "#FD8D3C" ,"#FC4E2A", 
    "#CCECE6", "#99D8C9", "#66C2A4","#41AE76", "#238B45",
    brewer.pal(9, "Greys")[4:8]
  ),order = rev(c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC","Endo","Acinar","Pericyte","EMT_like","Other"))
  )+theme_bw()+
    theme(plot.background = element_blank(),panel.background = element_blank(),panel.grid = element_blank())
  # pdf(file=paste0("A1_scRNA-seq_UMAP_",Sys.Date(),".pdf"),width=7,height = 5)
  # print(p1)
  # dev.off()
  
  Idents(TotalTissue.harmony) = TotalTissue.harmony$David_Tuevson_Clusters
  p2=DimPlot(TotalTissue.harmony,label = TRUE,raster=T,cols = c(
    "#F768A1","#9EBCDA",#"#8C96C6", "#8C6BB1", 
    "#FED976","#FC4E2A", #"#FEB24C", #"#FD8D3C" ,
    "#CCECE6", "#99D8C9","#41AE76", "#238B45",# "#66C2A4",
    brewer.pal(9, "Greys")[4:7]
  ))+theme_bw()+
    theme(plot.background = element_blank(),panel.background = element_blank(),panel.grid = element_blank())
  
  plot = ggarrange(p1,p2,ncol=2)
  pdf(file=paste0("A_scRNA-seq_UMAP_",Sys.Date(),".pdf"),width=14,height = 5)
  print(plot)
  dev.off()
  
  # df = FindAllMarkers(subcell.Seurat,subset.ident = c(8,9,10,13))
  write.csv(file="umap.csv",x = TotalTissue.harmony@meta.data)
  
}



#panel b, gene numbers, scRNA-seq
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_2022-09-14.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  
  ##> scRNA-seq
  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
  scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
  
  scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
  scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(Obj.list$meta$CellType))
  
  ##> plot
  TotalTissue.harmony@meta.data = scRNA_meta_filter.df
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels_V2
  
  pdf(file=paste0("B_scRNA-seq_GeneNumber_",Sys.Date(),".pdf"),width=7,height = 5)
  Seurat::VlnPlot(TotalTissue.harmony,features = "nFeature_RNA") + theme_bw() +
    scale_fill_manual(values = levels(Obj.list$meta$CellType_Color)) +
    theme(plot.background = element_blank(),panel.background = element_blank(),strip.background = element_blank(),
          legend.background = element_blank(),
          text = element_text(size=12),panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour = levels(Obj.list$meta$CellType_Color))
          
          )+
    labs(title="scRNA-seq")+xlab("")+ylab("Number of Genes")+
    RotatedAxis()
  dev.off()
}

#panel c, cell number 
{
  #V1, samples
  {
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_2022-09-14.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  
  ##> scRNA-seq
  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  scRNA_meta_filter.df = scRNA_meta_filter.df %>% select(sample,Main_subset_labels_V2)
  scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
  scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
  
  scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
  scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(Obj.list$meta$CellType))
  
  scRNA_number.df = scRNA_meta_filter.df %>% group_by(sample,Main_subset_labels_V2) %>% summarize(value=n())
  scRNA_number.df$omic = "scRNA"
  names(scRNA_number.df) = c("Sample","CellType","Value","Omic")
  ##> proteomics
  FACS_percentage.df = read.delim(file="FACS_percentage.txt",header = T) %>% t() %>% as.data.frame()
  colnames(FACS_percentage.df) = paste0("KPC",1:5)
  
  FACS_percentage_gather.df = tidyr::gather(FACS_percentage.df,key="Sample",value="Value")
  FACS_percentage_gather.df$CellType = rep(row.names(FACS_percentage.df),5)
  
  FACS_percentage_gather.df$Omic = "Proteomics"
  
  ##> combine
  plot.data =rbind(scRNA_number.df,FACS_percentage_gather.df)
  plot.data$CellType = factor(plot.data$CellType,levels=levels(Obj.list$meta$CellType))
  
  pdf(file=paste0("A_Comparision_RNAvsProteomics_",Sys.Date(),".pdf"),width=7,height = 5)
  ggplot(plot.data, aes(fill=CellType, y=Value, x=Sample)) + 
    geom_bar(position="fill", stat="identity") + 
    theme_classic() + facet_wrap(.~Omic) + xlab("")+ylab("") +
    scale_fill_manual(values = levels(Obj.list$meta$CellType_Color)) +
    theme(plot.background = element_blank(),panel.background = element_blank(),strip.background = element_blank(),
          legend.background = element_blank(),text = element_text(size=12),panel.grid = element_blank())#+
  # scale_y_continuous(expand = c(0,0))
  dev.off()
  
  write.csv(plot.data,file=paste0("A_Comparision_RNAvsProteomics_",Sys.Date(),".csv"))
  }
  
  #V2, combined, not used
  {
    rm(list=ls())
    TotalTissue.harmony = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_2022-09-14.rda")
    Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
    
    ##> scRNA-seq
    scRNA_meta.df = TotalTissue.harmony@meta.data
    scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
    scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
    scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
    scRNA_meta_filter.df = scRNA_meta_filter.df %>% select(sample,Main_subset_labels_V2)
    scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
    scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
    
    scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
    scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(Obj.list$meta$CellType))
    
    scRNA_number.df = scRNA_meta_filter.df %>% group_by(sample,Main_subset_labels_V2) %>% summarize(value=n())
    # scRNA_number.df$omic = "scRNA"
    names(scRNA_number.df) = c("Sample","CellType","Value")
    
    ## data transform
    scRNA_number_short.df = scRNA_number.df %>% spread(key="Sample",value="Value")
    
    col_sum = as.numeric(apply(scRNA_number_short.df[,2:5],2,sum))
    scRNA_percent_short.df=apply(scRNA_number_short.df[,2:5],1,FUN = function(x){
      x/col_sum*(100)
    }) %>% t() %>% as.data.frame()
    
    scRNA_percent_short.df$mean = rowMeans(scRNA_percent_short.df)
    scRNA_percent_short.df$sd = apply(scRNA_percent_short.df[,1:4],1,sd)*(1)
    
    # row.names(scRNA_percent_short.df) = scRNA_number_short.df$CellType
    scRNA_percent_short.df$CellType = scRNA_number_short.df$CellType
    
    scRNA_percent.df = scRNA_percent_short.df %>% #select(CellType,mean,sd) %>% 
      mutate(Omic = "scRNA_seq")
    write.csv(x =scRNA_percent.df, file=paste0("scRNA_percent_",Sys.Date(),".csv"))
    # scRNA_percent_long.df = scRNA_percent_short.df %>% 
    #   gather(key="Sample",value="Value",mean:sd)# %>% 
      # arrange(factor(CellType, levels = levels(Obj.list$meta$CellType))) #%>% 
      # mutate(x=factor(x, levels=unique(x)))
    # scRNA_percent_long.df$Omic = "scRNA_seq"
    
    ## not use
    ##> proteomics
    FACS_percentage.df = read.delim(file="FACS_percentage.txt",header = T) %>% t() %>% as.data.frame()
    colnames(FACS_percentage.df) = paste0("KPC",1:5)
    FACS_percentage.df$mean = rowMeans(FACS_percentage.df)
    FACS_percentage.df$sd = apply(FACS_percentage.df[,1:5],1,sd)
    
    # FACS_percentage_gather.df = tidyr::gather(FACS_percentage.df,key="Sample",value="Value")
    # FACS_percentage_gather.df$CellType = rep(row.names(FACS_percentage.df),5)
    FACS_percentage.df = FACS_percentage.df %>% select(mean,sd) %>% 
      mutate(Omic = "FACS")
    FACS_percentage.df$CellType = row.names(FACS_percentage.df)
    
    ##> combine
    plot.data =rbind(scRNA_percent.df,FACS_percentage.df)
    ggplot(plot.data, aes(fill=Omic, y=mean, x=CellType)) + 
      geom_bar(position="stack", stat="identity")#+
  
  }
}
#panel d, comparison of markers
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2022-09-14.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  
  ##> scRNA-seq
  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
  scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
  
  scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
  scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(Obj.list$meta$CellType))
  
  ##> plot
  scRNA_meta.df$Main_subset_labels_V2=factor(scRNA_meta.df$Main_subset_labels_V2,levels=(
    c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC"#,
      # "Endo","Acinar","Pericyte","Other"
      )))
  # levels(Obj.list$meta$CellType_Color)
  # levels(Obj.list$meta$CellType)
  # ))
  TotalTissue.harmony@meta.data = scRNA_meta.df
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels
  deg_genes = c(
    "Epcam", "Cdh1", "Tspan8", "Spp1", "Agr2", "Msln", "Krt19", 
    "Ptgis", "Tagln", "Col15a1", "Sparc", "Col8a1", "Fap", "Col14a1", "Tnc", "Thbs2", "Dcn", "Eng", "Pdgfra", "Serpine2", "Vim", 
    "Cd4", "Cd8a", "Cd8b1", "Foxp3", "Cd79a", "Ms4a1", "Cd79b", "Cd19", 
    "Itgam", "Ear2","Cd14", "Cd177", "S100a8", 
      "Adgre1", "Itgae","Cd163", "H2-Eb1", "H2-Ab1", "H2-Aa", "Cd74", "Itgax"
    
  )
  set.seed(1000)
  pdf(file=paste0("D_scRNA-seq_Pheatmap_",Sys.Date(),".pdf"),width=5,height = 7)
  DoHeatmap(subset(TotalTissue.harmony, downsample = 100,idents = c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC")),
            features = deg_genes, size = 3,draw.lines = F,raster = F
            )+scale_fill_gradientn(colors=c("blue","white","red"))
  dev.off()
}

#panel d, comparison of markers, V2
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_2022-09-14.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  
  ##> scRNA-seq
  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
  scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
  
  scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
  scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(Obj.list$meta$CellType))
  
  ##> plot
  scRNA_meta.df$Main_subset_labels_V2=factor(scRNA_meta.df$Main_subset_labels_V2,levels=(
    c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC"#,
      # "Endo","Acinar","Pericyte","Other"
    )))
  # levels(Obj.list$meta$CellType_Color)
  # levels(Obj.list$meta$CellType)
  # ))
  TotalTissue.harmony@meta.data = scRNA_meta.df
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels
  # deg_genes = c(
  #   "Epcam", "Cdh1", "Tspan8", "Spp1", "Agr2", "Msln", "Krt19", 
  #   "Ptgis", "Tagln", "Col15a1", "Sparc", "Col8a1", "Fap", "Col14a1", "Tnc", "Thbs2", "Dcn", "Eng", "Pdgfra", "Serpine2", "Vim", 
  #   "Cd4", "Cd8a", "Cd8b1", "Foxp3", "Cd79a", "Ms4a1", "Cd79b", "Cd19", 
  #   "Itgam", "Ear2","Cd14", "Cd177", "S100a8", 
  #   "Adgre1", "Itgae","Cd163", "H2-Eb1", "H2-Ab1", "H2-Aa", "Cd74", "Itgax"
  #   
  # )
  deg_genes = c("Epcam", "Msln", "Krt19","Pdpn","Ly6c1","Tagln","Il2ra","Foxp3", "Cd19",
                "Itgam","Cd74", "Itgax")
  set.seed(1000)
  pdf(file=paste0("D_scRNA-seq_Pheatmap_",Sys.Date(),".pdf"),width=5,height = 6)
  DoHeatmap(subset(TotalTissue.harmony, downsample = 100,idents = c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC")),
            features = deg_genes, size = 3,draw.lines = F,raster = F
  )+scale_fill_gradientn(colors=c("blue","white","red"))
  dev.off()
}

#panel d, comparison of markers, V3
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2022-11-29.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  
  ##> scRNA-seq
  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
  scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
  
  scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
  scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(Obj.list$meta$CellType))
  
  ##> plot
  scRNA_meta.df$Main_subset_labels_V2=factor(scRNA_meta.df$Main_subset_labels_V2,levels=(
    c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC"#,
      # "Endo","Acinar","Pericyte","Other"
    )))
  # levels(Obj.list$meta$CellType_Color)
  # levels(Obj.list$meta$CellType)
  # ))
  TotalTissue.harmony@meta.data = scRNA_meta.df
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels
  deg_genes = c(
    "Epcam", "Cdh1", "Tspan8", "Spp1", "Agr2", "Msln", "Krt19",
    "Ptgis", "Tagln", "Col15a1", "Sparc", "Col8a1", "Fap", "Col14a1", "Tnc", "Thbs2", "Dcn", "Eng", "Pdgfra", "Serpine2", "Vim",
    "Cd4", "Cd8a", "Cd8b1", "Foxp3", "Cd79a", "Ms4a1", "Cd79b", "Cd19",
    "Itgam", "Ear2","Cd14", "Cd177", "S100a8",
    "Adgre1", "Itgae","Cd163", "H2-Eb1", "H2-Ab1", "H2-Aa", "Cd74", "Itgax"

  )
  deg_genes2 = c(
    "Epcam", "Cdh1", "Tspan8", "Spp1", "Agr2", "Msln", "Krt19",
    "Ptgis", "Tagln", "Col15a1", "Sparc", "Col8a1", "Fap", "Col14a1", "Tnc", "Thbs2", "Dcn", "Eng", "Pdgfra", "Serpine2", "Vim",
    "Cd4", "Cd8a", "Cd8b", "Foxp3", "Cd79a", "Ms4a1", "Cd79b", "Cd19",
    "Itgam", "Ear2","Cd14", "Cd177", "S100a8",
    "Adgre1", "Itgae","Cd163", "H2-Eb1", "H2-Ab1", "H2-Aa", "Cd74", "Itgax"
    
  )
  # deg_genes = c("Epcam", "Msln", "Krt19","Pdpn","Ly6c1","Tagln","Il2ra","Foxp3", "Cd19",
  #               "Itgam","Cd74", "Itgax")
  ##>RNA
  RNA_average.df = AverageExpression(TotalTissue.harmony,features = deg_genes)$RNA %>% as.data.frame()
  
  ##>Protein
  # copy_num.df = read.delim(file="../../resource/intensity_copynum.txt",header = T)
  # copy_num.df = copy_num.df[,c(141,142,c(72:140))]
  # row.names(copy_num.df) = copy_num.df$Protein.ID
  # copy_num.df$Protein.ID = NULL
  # copy_num.df$Gene = NULL
  # names(copy_num.df) = gsub(names(copy_num.df),pattern = "^Copy.number.|_[0-9]+.Total.Intensity",replacement = "")
  # copy_num.df = copy_num.df[,Obj.list$meta$RawSampleID]
  # names(copy_num.df) = Obj.list$meta$SampleID
  lfq.df = Obj.list$filter
  
  meta.df = Obj.list$meta
  lfq_mean.df = data.frame(row.names = row.names(lfq.df))
  for(i in levels(meta.df$CellType)[-2]){
    temp_index = which(meta.df$CellType==i)  
    lfq_mean.df[,i] = rowMeans(lfq.df[,temp_index],
                               na.rm = TRUE)
  }
  
  anno.df = Obj.list$annotation[Obj.list$annotation$genename %in% deg_genes2,]
  anno.df = anno.df[c(1:39,41,42,44),]
  lfq_mean.df = lfq_mean.df[anno.df$pid,]
  row.names(lfq_mean.df) = anno.df$genename
  lfq_mean.df = lfq_mean.df[deg_genes2,]
  
  #combine
  RNA_average.df$order = seq(from=1,to=83,by=2)
  
  lfq_mean.df$order = seq(from=2,to=84,by=2)
  RNA_average.df = RNA_average.df %>% select(names(lfq_mean.df))
  
  plot.df = rbind(RNA_average.df,lfq_mean.df)
  plot.df = plot.df[order(plot.df$order),]
  plot.df$order = NULL
  
  plot_anno.df = data.frame(group=factor(rep(c("RNA","Protein"),42),levels = c("RNA","Protein")))
  row.names(plot_anno.df) = row.names(plot.df)
  # row.names(plot.df) = gsub(row.names(plot.df),pattern = "-",replacement = "_")
  pdf(file=paste0("Fig4_D_RNA-Protein-DEG_Pheatmap_",Sys.Date(),".pdf"),width=5,height = 8)
  pheatmap::pheatmap(plot.df,cellwidth = 18,cellheight = 6,
                     scale = "row",border_color = "white",
                     cluster_rows=FALSE,cluster_cols = FALSE,angle_col=45,
                     annotation_names_row = F,annotation_names_col = F,legend = T,
                     fontsize = 5,display_numbers = FALSE,
                     #silent = FALSE,
                     annotation_row  = plot_anno.df,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                     # color = colorRampPalette(c("blue","white","red"))(100), 
                     width = 9,height = 5)
  dev.off()
  
  mat = as.matrix(plot.df)
  mat = apply(mat,1,scale) %>% t() %>% as.data.frame()
  names(mat) = names(plot.df)
  write.csv(file="scRNA-seq_pheatmap.csv",x = mat)
  # copy_num.df$copy_num_mean = apply(copy_num.df[,3:71],1,sum)
  # copy_num.df = copy_num.df %>% select(Protein.ID,Gene,copy_num_mean) %>% filter(copy_num_mean!=0)
  
  # set.seed(1000)
  # pdf(file=paste0("D_scRNA-seq_Pheatmap_",Sys.Date(),".pdf"),width=5,height = 6)
  # DoHeatmap(subset(TotalTissue.harmony, downsample = 100,idents = c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC")),
  #           features = deg_genes, size = 3,draw.lines = F,raster = F
  # )+scale_fill_gradientn(colors=c("blue","white","red"))
  # dev.off()
}
#panel supp, used markers
{
  rm(list=ls())
  used_markers=c(
  "Epcam", "Krt19", "Krt8",
  "Cdh5", "Plvap", "Vwf", #Endo,16
  "Ctrb1", "Prss2","Try5", # Acinar,17
  "Rgs5", "Acta2", "Pdgfrb", #Pericyte,18
  "Lum", "Dcn", "Col1a1", "Col3a1", #Fib,11
  "Cd4", "Cd8a","Foxp3",
  "Cd79a","Cd19",
  "Ear2","Lyz1",
  "S100a8","S100a9","G0s2",
  #"Cd14","Csf1r","Itgam",
  # "Fcgr1","Cx3cr1","Treml4", # Mono
  "Adgre1","Saa3","Spp1","C1qa","C1qb","C1qc",#"Marco", #1,4
  #"Fscn1", "Ccr7",
  #"Ccl22","Cd274","Ccl5",
  "H2-Ab1","H2-Eb1","Cd74"
  )

  TotalTissue.harmony = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_2022-09-14.rda")
  # Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  # TotalTissue.harmony$Main_cell_labels = factor(TotalTissue.harmony$Main_cell_labels,levels = c(
  #   "PCC","Acinar","Endo","Pericyte","CAF","T","B","MYE","NEU",
  # ))
  Idents(TotalTissue.harmony)=TotalTissue.harmony$Main_subset_labels
  
  pdf(file=paste0("supp_A_scRNA-seq_UsedMarker_",Sys.Date(),".pdf"),width=11,height = 7)
  DotPlot(TotalTissue.harmony,features = used_markers)+RotatedAxis()+xlab("")+ylab("")
  dev.off()
 
}

#panel, FACS used markers; Vlnplot, correlation plot
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2022-09-14.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  
  ##> scRNA-seq
  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  
  scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels = c(
    "PCC","iCAF","myCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC"
  ))
  TotalTissue.harmony@meta.data = scRNA_meta_filter.df
  TotalTissue.harmony.filter = subset(TotalTissue.harmony,
                                      idents=unique(scRNA_meta_filter.df$Main_subset_labels_V2))
  
  Idents(TotalTissue.harmony.filter) = TotalTissue.harmony.filter$Main_subset_labels_V2
  
  saveRDS(TotalTissue.harmony.filter,file = paste0("GSE129455_All_Viable_with_subset_filter_",Sys.Date(),".rda"))
  
  
  ##> Vlnplot
  p1=VlnPlot(TotalTissue.harmony.filter,features = c(
    "Epcam","Msln","Pdpn",
    "Ly6c1","Cd74",
    "Ptprc","Cd3e","Cd4","Cd8a","Il2ra","Cd19",
    "Itgam","Ly6g","Adgre1","Itgax"),ncol = 3,combine = TRUE,pt.size = 0,same.y.lims = FALSE,
          cols = levels(Obj.list$meta$CellType_Color)[-2]) +RotatedAxis()#+ theme_bw()+
          # theme(plot.background = element_blank(),panel.grid = element_blank(),
          #       legend.position="none",plot.title = element_text(hjust = 0.5,size = 14))+
          # labs(y="",fill="",x="")
  pdf(file=paste0("H_scRNA_FACS_marker_",Sys.Date(),".pdf"),width=10,height=15)
  print(p1)
  dev.off()
  ##> corr, RNA
  av.exp <- AverageExpression(TotalTissue.harmony.filter)$RNA
  cor.exp <- as.data.frame(cor(av.exp))
  
  library(corrplot)
  library(GGally)
  
  
  p1=ggcorr(cor_matrix =cor(av.exp)[],data=NULL, midpoint=0.5,label = FALSE,
         label_round=2,label_size = 4,size=5,low="steelblue",mid="white",high="red",
         legend.position = NULL,limits=c(0,1),label_alpha = TRUE)+theme(legend.position = "none")#+coord_flip()
  write.csv(cor(av.exp)[],file="rna_cor.csv")
  ##> protein
  meta.df = Obj.list$meta
  lfq_mean.df = data.frame(row.names = row.names(Obj.list$filter))
  for(i in levels(meta.df$CellType)[-2]){
    temp_index = which(meta.df$CellType==i)  
    lfq_mean.df[,i] = rowMeans(Obj.list$filter[,temp_index],
                                            na.rm = TRUE)
  }#darkgrey","white","darkgoldenrod1","brown3"
  p2=ggcorr(cor_matrix =cor(lfq_mean.df)[],data=NULL, midpoint =0.5,label = FALSE,
         label_round=2,label_size = 4,size=5,low="steelblue",mid="white",high="red",
         legend.position = NULL,limits=c(0,1),label_alpha = TRUE)+theme(legend.position = "none")+coord_flip()
  write.csv(cor(lfq_mean.df)[],file="protein_cor.csv")
  pdf(file=paste0("F_corr_protein_",Sys.Date(),".pdf"))
  print(p2)
  dev.off()
  pdf(file=paste0("F_corr_rna_",Sys.Date(),".pdf"))
  print(p1)
  dev.off()
  
  
  
}

#panel, FACS used markers; Vlnplot, V2
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_2022-09-14.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  
  ##> scRNA-seq
  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  
  scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels = c(
    "PCC","iCAF","myCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC"
  ))
  TotalTissue.harmony@meta.data = scRNA_meta_filter.df
  TotalTissue.harmony.filter = subset(TotalTissue.harmony,
                                      idents=unique(scRNA_meta_filter.df$Main_subset_labels_V2))
  
  Idents(TotalTissue.harmony.filter) = TotalTissue.harmony.filter$Main_subset_labels_V2
  
  saveRDS(TotalTissue.harmony.filter,file = paste0("scRNA/GSE129455_All_Viable_with_subset_filter_",Sys.Date(),".rda"))
  
  
  ##> Vlnplot
  p1=VlnPlot(TotalTissue.harmony.filter,features = c("Epcam", "Msln", "Krt19","Pdpn","Ly6c1","Tagln","Il2ra","Foxp3", "Cd19",
                                                                   "Itgam","Cd74", "Itgax"),
             ncol = 3,combine = TRUE,pt.size = 0,same.y.lims = FALSE,
    cols = levels(Obj.list$meta$CellType_Color)[-2]) +RotatedAxis()#+ theme_bw()+
  # theme(plot.background = element_blank(),panel.grid = element_blank(),
  #       legend.position="none",plot.title = element_text(hjust = 0.5,size = 14))+
  # labs(y="",fill="",x="")
  pdf(file=paste0("H_scRNA_FACS_marker_",Sys.Date(),".pdf"),width=10,height=12)
  print(p1)
  dev.off()
  ##> corr, RNA
  av.exp <- AverageExpression(TotalTissue.harmony.filter)$RNA
  cor.exp <- as.data.frame(cor(av.exp))
  
  library(corrplot)
  library(GGally)
  
  
  p1=ggcorr(cor_matrix =cor(av.exp)[],data=NULL, midpoint=0.5,label = FALSE,
            label_round=2,label_size = 4,size=5,low="steelblue",mid="white",high="red",
            legend.position = NULL,limits=c(0,1),label_alpha = TRUE)+theme(legend.position = "none")#+coord_flip()
  
  ##> protein
  meta.df = Obj.list$meta
  lfq_mean.df = data.frame(row.names = row.names(Obj.list$filter))
  for(i in levels(meta.df$CellType)[-2]){
    temp_index = which(meta.df$CellType==i)  
    lfq_mean.df[,i] = rowMeans(Obj.list$filter[,temp_index],
                               na.rm = TRUE)
  }#darkgrey","white","darkgoldenrod1","brown3"
  ggcorr(cor_matrix =cor(lfq_mean.df)[],data=NULL, midpoint =0.5,label = FALSE,
         label_round=2,label_size = 4,size=5,low="steelblue",mid="white",high="red",
         legend.position = NULL,limits=c(0,1),label_alpha = TRUE)+theme(legend.position = "none")+coord_flip()
  
  pdf(file=paste0("F_corr_protein_",Sys.Date(),".pdf"))
  print(p2)
  dev.off()
  pdf(file=paste0("F_corr_rna_",Sys.Date(),".pdf"))
  print(p1)
  dev.off()
  
  print(p2)
  dev.off()
  
  
}

#panel, selected markers;bar plot
{
  rm(list=ls())
  deg_genes = c(
    "Epcam","Msln","Tagln",
    "Foxp3","Cd19","Cd74",
    "Itgae","Itgax"
    
  )
  
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_filter_2022-09-15.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  meta.df = Obj.list$meta
  
  anno.df = Obj.list$annotation[Obj.list$annotation$genename %in% deg_genes,]
  row.names(anno.df) = anno.df$genename
  anno.df = anno.df[deg_genes,]#sort
  ##> scRNA-seq
  
  ##> RNA
  log2_RNA_average = AverageExpression(TotalTissue.harmony,features = deg_genes,slot = "counts")$RNA %>% as.data.frame()
  ##> Protein
  log2_lfq_average.df = data.frame(row.names = row.names(Obj.list$filter))
  for(i in levels(meta.df$CellType)[-2]){
    temp_index = which(meta.df$CellType==i)  
    log2_lfq_average.df[,i] = rowMeans(log2(Obj.list$filter[,temp_index]+1),
                               na.rm = TRUE)
  }#darkgrey","white","darkgoldenrod1","brown3"
  log2_lfq_average.df = log2_lfq_average.df[anno.df$pid,]
  row.names(log2_lfq_average.df) = anno.df$genename
  ##> Combine
  write.csv(log2_RNA_average,file="log2_RNA_average.csv")
  write.csv(log2_lfq_average.df,file="log2_lfq_average.csv")
  
  # plot.data1 = log2_RNA_average %>% as.data.frame()
  plot.data1 = log2_RNA_average %>% apply(1,FUN = function(x){scale(x,center = F)}) %>% t() %>% as.data.frame()
  names(plot.data1) = names(log2_RNA_average)
  plot.data1$type = "RNA"
  plot.data1$genename = row.names(plot.data1)
  plot.data1.gather = plot.data1 %>% gather(key="CellType",value="log2_exp",-type,-genename)
  plot.data1.gather$log2_exp = -plot.data1.gather$log2_exp
  
  plot.data2 = log2_lfq_average.df %>% apply(1,FUN = function(x){scale(x,center = F)}) %>% t() %>% as.data.frame()
  names(plot.data2) = names(log2_lfq_average.df)
  plot.data2$type = "Protein"
  plot.data2$genename = row.names(plot.data2)
  plot.data2.gather = plot.data2 %>% gather(key="CellType",value="log2_exp",-type,-genename)
  
  plot.data = rbind(plot.data1.gather,plot.data2.gather)
  plot.data$CellType = factor(plot.data$CellType,levels = levels(Obj.list$meta$CellType)[-2])
  
  write.csv(plot.data,file="relative_exp.csv")#output
  
  plot.list = list()
  for(i in unique(plot.data$genename)){
    temp_plot.data = plot.data %>% filter(genename==i)
    p1=ggplot(data=temp_plot.data) + geom_col(aes(x=CellType,y=log2_exp,fill=type))    +
      scale_y_continuous(breaks = seq(from = -10, to = 10, by = 1), labels = c(seq(10,0, -1), seq(1, 10, 1))) + 
      theme_bw()+
      theme(text = element_text(size = 10),legend.title=element_blank(),title = element_text(hjust = 0.5))+
      labs(x="",y="Reletive abundance",title="")+
      # guides(fill="none")+
      RotatedAxis() #+ coord_flip() 
    plot.list[[i]] = p1
  }
  library(ggpubr)
  p1=ggarrange(plotlist = plot.list,ncol=3,common.legend=T,legend="right",nrow = 3,labels = unique(plot.data$genename))
  
  pdf(file="Fig4_deg_genes_RNA-Proteins_2.pdf",width = 7,height = 7)
  print(p1)
  dev.off()
  
}

#panel, GO enrichment analysis.
{
  rm(list=ls())
  TotalTissue_filter.seurat = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_filter_2022-09-15.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  deg.df = FindAllMarkers(TotalTissue_filter.seurat,logfc.threshold = 1)
  deg_cutoff.df = deg.df %>% filter(p_val_adj<=0.05 ) %>% filter(avg_log2FC>=1)
  deg_cutoff.df$CellLineage = "Myeloid"
  deg_cutoff.df[deg_cutoff.df$cluster %in% "PCC","CellLineage"]="Cancer Cell"
  deg_cutoff.df[deg_cutoff.df$cluster %in% c("iCAF","myCAF","apCAF"),"CellLineage"]="CAF"
  deg_cutoff.df[deg_cutoff.df$cluster %in% c("T4","T8","Treg","B"),"CellLineage"]="Lymphoid Cell"
  
  ##> GO analysis
  library(DOSE)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  genename1_tran.df <- bitr(deg_cutoff.df$gene, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  
  deg_cutoff.df = deg_cutoff.df %>% merge(genename1_tran.df,by.x="gene",by.y="SYMBOL",all.x=FALSE,all.y=TRUE) ##883
  
  
  genename_bg_tran.df <- bitr(row.names(TotalTissue_filter.seurat[["RNA"]]), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  
  deg_cutoff.df$CellLineage = factor(deg_cutoff.df$CellLineage,levels=c(
    c("Cancer Cell","CAF","Lymphoid Cell","Myeloid")
  ))
  
  ##with CellLineage
  formula_res <- compareCluster(ENTREZID~CellLineage, data=deg_cutoff.df, fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,
                                universe = genename_bg_tran.df$ENTREZID,readable=TRUE
  )
  
  formula_res_cutoff = formula_res
  formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<=0.05,]
  formula_res_cutoff@compareClusterResult$CellLineage = factor(formula_res_cutoff@compareClusterResult$CellLineage,
                                                               levels=c(
    c("Cancer Cell","CAF","Lymphoid Cell","Myeloid")
  ))
  p1=dotplot(formula_res_cutoff,x="CellLineage",color="pvalue",showCategory = 3,label_format=35,font.size=12) 
  
  p2 = p1
  p2$theme$axis.text.x$angle=45
  p2$theme$axis.text.x$hjust = 0.5
  p2$theme$axis.text.x$vjust = 0.5
  p2$theme$axis.text.x$colour= levels(Obj.list$meta$CellLineage_Color)
  p2$labels$x=NULL
  pdf(file=paste0("Fig4_scRNA-seq_GOBP_CellLineage_",Sys.Date(),".pdf"),
      width=6,height=6)
  p2
  dev.off()
  write.csv(formula_res_cutoff@compareClusterResult,file=paste0("Fig4_scRNA-seq_GOBP_CellLineage_",Sys.Date(),".csv"))
  
  
  
}


#panel, CAF de-convlution
{
  rm(list=ls())
  TotalTissue_filter.seurat = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_filter_2022-09-15.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  deg.df = FindAllMarkers(TotalTissue_filter.seurat,logfc.threshold = 1)
  # write.csv(x=deg.df,file=paste0("findAllMarkers_scRNA_",Sys.Date(),".csv"))
  
  deg_cutoff.df = deg.df %>% filter(p_val_adj<=0.05 ) %>% filter(avg_log2FC>=1)
  deg_cutoff.df$CellLineage = "Myeloid"
  deg_cutoff.df[deg_cutoff.df$cluster %in% "PCC","CellLineage"]="Cancer Cell"
  deg_cutoff.df[deg_cutoff.df$cluster %in% c("iCAF","myCAF","apCAF"),"CellLineage"]="CAF"
  deg_cutoff.df[deg_cutoff.df$cluster %in% c("T4","T8","Treg","B"),"CellLineage"]="Lymphoid Cell"
  
  ##> CAF
  CAF_deg_cutoff.df = deg_cutoff.df %>% filter(CellLineage == "CAF") %>% unique()
  
  ## CAF barcode and average sig matrix
  CAF_TotalTissue.seurat = subset(TotalTissue_filter.seurat,
                                  features = unique(CAF_deg_cutoff.df$gene),idents = c("iCAF","myCAF","apCAF"))
  CAF_TotalTissue.seurat$Main_subset_labels_V2 = factor(CAF_TotalTissue.seurat$Main_subset_labels_V2,levels=(
    c("iCAF","myCAF","apCAF")))
  CAF_deg.df = AverageExpression(CAF_TotalTissue.seurat,slot="counts")$RNA %>% data.frame()
  CAF_deg.df$`Gene symbol` = row.names(CAF_deg.df)
  write.csv(CAF_deg.df,file=paste0("CAF_subtype_markers_",Sys.Date(),".csv"),row.names = F)
}

#panel, CAF de-convlution V2
{
  rm(list=ls())
  TotalTissue_filter.seurat = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_filter_2022-09-15.rda")
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  deg.df = FindAllMarkers(TotalTissue_filter.seurat,logfc.threshold = 1)
  # write.csv(x=deg.df,file=paste0("findAllMarkers_scRNA_",Sys.Date(),".csv"))
  
  deg_cutoff.df = deg.df %>% filter(p_val_adj<=0.05 ) %>% filter(avg_log2FC>=1)
  deg_cutoff.df$CellLineage = "Myeloid"
  deg_cutoff.df[deg_cutoff.df$cluster %in% "PCC","CellLineage"]="Cancer Cell"
  deg_cutoff.df[deg_cutoff.df$cluster %in% c("iCAF","myCAF","apCAF"),"CellLineage"]="CAF"
  deg_cutoff.df[deg_cutoff.df$cluster %in% c("T4","T8","Treg","B"),"CellLineage"]="Lymphoid Cell"
  
  ##> CAF
  CAF_deg_cutoff.df = deg_cutoff.df %>% filter(CellLineage == "CAF") %>% unique()
  
  CAF_deg = names(table(CAF_deg_cutoff.df$gene)[which(table(CAF_deg_cutoff.df$gene)==1)])
  ## CAF barcode and average sig matrix
  CAF_TotalTissue.seurat = subset(TotalTissue_filter.seurat,
                                  features = CAF_deg,idents = c("iCAF","myCAF","apCAF"))
  CAF_TotalTissue.seurat$Main_subset_labels_V2 = factor(CAF_TotalTissue.seurat$Main_subset_labels_V2,levels=(
    c("iCAF","myCAF","apCAF")))
  CAF_deg.df = AverageExpression(CAF_TotalTissue.seurat,slot="counts")$RNA %>% data.frame()
  CAF_deg.df$`Gene symbol` = row.names(CAF_deg.df)
  write.csv(CAF_deg.df,file=paste0("CAF_subtype_markers_V2_",Sys.Date(),".csv"),row.names = F)
}


##Cor V1,all vs all
{
rm(list=ls())
copy_num.df = read.delim(file="intensity_copynum.txt",header = T)
copy_num.df = copy_num.df[,c(141,142,c(72:140))]
row.names(copy_num.df) = copy_num.df$Protein.ID
copy_num.df$copy_num_mean = apply(copy_num.df[,3:71],1,sum)
copy_num.df = copy_num.df %>% select(Protein.ID,Gene,copy_num_mean) %>% filter(copy_num_mean!=0)

TotalTissue_filter.seurat = readRDS(file = "scRNA/GSE129455_All_Viable_with_subset_filter_2022-09-15.rda")
df =exp(TotalTissue_filter.seurat[["RNA"]]@data)-1
average.df = apply(df,1,sum)
# average.df = AverageExpression(TotalTissue_filter.seurat,slot = "counts")$RNA
scRNA_mean.df = average.df
scRNA_mean.df = scRNA_mean.df %>% as.data.frame()
# scRNA_mean.df$relative_exp = scRNA_mean.df$. / max(scRNA_mean.df$.)
scRNA_mean.df$genename = row.names(scRNA_mean.df)
names(scRNA_mean.df) = c("exp","genename")
scRNA_mean.df = scRNA_mean.df %>% select(genename,exp) %>% filter(exp!=0)

Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
anno.df = Obj.list$annotation %>% filter(genename!="")

comparison.df = anno.df %>% merge(copy_num.df,by.x="pid",by.y="Protein.ID",all.x=F,all.y=F)
comparison.df = comparison.df %>% merge(scRNA_mean.df,by="genename",all.x=F,all.y=F)
comparison.df$relative_copynum = comparison.df$copy_num_mean/max(comparison.df$copy_num_mean)
comparison.df$relative_exp = comparison.df$exp/max(comparison.df$exp)
df = comparison.df %>% arrange(desc(relative_copynum)) %>% 
  mutate(relative_copynum_rank=1:dim(comparison.df)[1]) %>% 
  arrange(desc(relative_exp)) %>% 
  mutate(relative_exp_rank=1:dim(comparison.df)[1])

df$copy_num_log2 = log2(df$copy_num_mean)
df$copy_num_zscore = (df$copy_num_log2-mean(df$copy_num_log2))/sd(df$copy_num_log2)

df$exp_log2 = log2(df$relative_exp)
df$exp_zscore = (df$exp_log2-mean(df$exp_log2))/sd(df$exp_log2)
# comparison.df$log2_copynum = log2(comparison.df$relative_copynum)
# comparison.df$log2_exp = log2(comparison.df$relative_exp)
# df = df %>% filter(copy_num_mean!=0 & exp!=0)

library(viridis)
p1=ggplot(data = df, mapping = aes(x =exp_zscore, y = copy_num_zscore)) +
  geom_bin2d(bins = 60) + # bins控制着图中每个箱子的大小
  scale_fill_viridis() +
  labs(tag = "Pearson.Corr=0.473") +
  theme_bw()+theme(plot.background = element_blank(),
                   panel.background = element_blank(),legend.background = element_blank(),panel.grid = element_blank()
  )
cor.test(df$copy_num_zscore,df$exp_zscore,method = "pearson",exact=TRUE)

p2=ggplot(data = df, mapping = aes(x =-log10(relative_exp_rank), y = -log10(relative_copynum_rank))) +
  geom_bin2d(bins = 60) + # bins控制着图中每个箱子的大小
  scale_fill_viridis() +
  labs(tag = "Pearson.Corr=0.500") +
  theme_bw()+theme(plot.background = element_blank(),
                   panel.background = element_blank(),legend.background = element_blank(),panel.grid = element_blank()
                    )
cor.test(df$relative_exp_rank,df$relative_copynum_rank,method = "pearson",exact=TRUE)

pdf(file=paste0("comparison_zscore_",Sys.Date(),".pdf"),width = 5,height = 3)
print(p1)
dev.off()

pdf(file=paste0("comparison_rank_",Sys.Date(),".pdf"),width = 5,height = 3)
print(p2)
dev.off()

write.csv(df,file=paste0("comparison_copynumvsexp_",Sys.Date(),".csv"))
temp.Seurat = readRDS(file="D:\\workspace\\01-PDAC\\05-FACS\\202209_Fig\\20220905\\Fig5\\scRNA\\T_cell\\38864\\10X_Treg_Scale_2022-09-23.rda")
}

##Cor V2,Cell vs Cell
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20220905.rds")
  # copy_num.df = read.delim(file="intensity_copynum.txt",header = T)
  # copy_num.df = copy_num.df[,c(141,142,c(72:140))]
  # row.names(copy_num.df) = copy_num.df$Protein.ID
  # copy_num.df$copy_num_mean = apply(copy_num.df[,3:71],1,sum)
  # copy_num.df = copy_num.df %>% select(Protein.ID,Gene,copy_num_mean) %>% filter(copy_num_mean!=0)
  
  copy_num.df = read.delim(file="../../resource/intensity_copynum.txt",header = T)
  copy_num.df = copy_num.df[,c(141,142,c(72:140))]
  row.names(copy_num.df) = copy_num.df$Protein.ID
  copy_num.df$Protein.ID = NULL
  copy_num.df$Gene = NULL
  names(copy_num.df) = gsub(names(copy_num.df),pattern = "^Copy.number.|_[0-9]+.Total.Intensity",replacement = "")
  copy_num.df = copy_num.df[,Obj.list$meta$RawSampleID]
  names(copy_num.df) = Obj.list$meta$SampleID
  # copy_num.df = log2(copy_num.df+1)
  
  
  meta.df = Obj.list$meta
  copynum_sum.df = data.frame(row.names = row.names(copy_num.df))
  for(i in levels(meta.df$CellType)[-2]){
    temp_index = which(meta.df$CellType==i)  
    copynum_sum.df[,i] = rowSums(copy_num.df[,temp_index],
                               na.rm = TRUE)
  }
  
  TotalTissue_filter.seurat = readRDS(file = "GSE129455_All_Viable_with_subset_filter_2022-09-15.rda")
  output_cor = c()
  for(i in levels(Obj.list$meta$CellType)[-2]){
  temp_RNA.df = subset(TotalTissue_filter.seurat,idents=c(i))
  temp_RNA.df =exp(temp_RNA.df[["RNA"]]@data)-1
  count_sum.df = apply(temp_RNA.df,1,sum) %>% as.data.frame()
  names(count_sum.df) = "RNA"
  count_sum.df$genename = row.names(count_sum.df)
  #merge
  temp_protein.df = copynum_sum.df %>% select(i)
  names(temp_protein.df) = "Protein"
  temp_protein.df$pid = row.names(copynum_sum.df)
  anno.df = Obj.list$annotation %>% filter(genename!="")
  comparison.df = anno.df %>% merge(temp_protein.df,by.x="pid",by.y="pid",all.x=F,all.y=F)
  comparison.df = comparison.df %>% merge(count_sum.df,by="genename",all.x=F,all.y=F)
  # comparison.df = comparison.df %>% filter("PCC"!=0 )
  comparison.df = comparison.df[-which(comparison.df$Protein==0 |comparison.df$RNA==0),]
  comparison.df$log2_Protein = log2(comparison.df$Protein)
  comparison.df$log2_RNA = log2(comparison.df$RNA)
  
  cor = cor.test(comparison.df$log2_Protein,comparison.df$log2_RNA,method = "pearson",exact=TRUE)
  output_cor = c(output_cor,paste0(i," ",as.numeric(cor$estimate)))
  }
  # 
  # df$copy_num_zscore = (df$copy_num_log2-mean(df$copy_num_log2))/sd(df$copy_num_log2)
  # 
  # df$exp_log2 = log2(df$relative_exp)
  # df$exp_zscore = (df$exp_log2-mean(df$exp_log2))/sd(df$exp_log2)
  # # comparison.df$log2_copynum = log2(comparison.df$relative_copynum)
  # # comparison.df$log2_exp = log2(comparison.df$relative_exp)
  # # df = df %>% filter(copy_num_mean!=0 & exp!=0)
  # 
  # library(viridis)
  # p1=ggplot(data = df, mapping = aes(x =exp_zscore, y = copy_num_zscore)) +
  #   geom_bin2d(bins = 60) + # bins控制着图中每个箱子的大小
  #   scale_fill_viridis() +
  #   labs(tag = "Pearson.Corr=0.473") +
  #   theme_bw()+theme(plot.background = element_blank(),
  #                    panel.background = element_blank(),legend.background = element_blank(),panel.grid = element_blank()
  #   )
  # cor.test(df$copy_num_zscore,df$exp_zscore,method = "pearson",exact=TRUE)
  # 
  # p2=ggplot(data = df, mapping = aes(x =-log10(relative_exp_rank), y = -log10(relative_copynum_rank))) +
  #   geom_bin2d(bins = 60) + # bins控制着图中每个箱子的大小
  #   scale_fill_viridis() +
  #   labs(tag = "Pearson.Corr=0.500") +
  #   theme_bw()+theme(plot.background = element_blank(),
  #                    panel.background = element_blank(),legend.background = element_blank(),panel.grid = element_blank()
  #   )
  # cor.test(df$relative_exp_rank,df$relative_copynum_rank,method = "pearson",exact=TRUE)
  # 
  # pdf(file=paste0("comparison_zscore_",Sys.Date(),".pdf"),width = 5,height = 3)
  # print(p1)
  # dev.off()
  # 
  # pdf(file=paste0("comparison_rank_",Sys.Date(),".pdf"),width = 5,height = 3)
  # print(p2)
  # dev.off()
  # 
  # write.csv(df,file=paste0("comparison_copynumvsexp_",Sys.Date(),".csv"))
  # temp.Seurat = readRDS(file="D:\\workspace\\01-PDAC\\05-FACS\\202209_Fig\\20220905\\Fig5\\scRNA\\T_cell\\38864\\10X_Treg_Scale_2022-09-23.rda")
}