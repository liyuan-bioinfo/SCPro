#'@author Yuan
#'@desc comparision with scRNA-seq
#'@version comparison and visualization of shared celltype using our ct_dataset_01 and public scRNA-seq dataset downloaded from GSE129455

library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyr)

library(corrplot)
library(GGally)
options(stringsAsFactors = F)
options(encoding = "UTF-8")

setwd("")

# Figure_13a, Comparision umap of shared cell types and origin cell types from source data's paper. 
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2024-03-27.rda") #Anyone can download processed scRNA-seq from https://doi.org/10.5281/zenodo.13981062 
  
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels2
  p1=DimPlot(TotalTissue.harmony,label = TRUE,raster = T,cols = c(
    "#F768A1","#9EBCDA","#8C96C6", "#8C6BB1", 
    "#FED976","#FEB24C", "#FD8D3C" ,"#FC4E2A", 
    "#CCECE6", "#99D8C9", "#66C2A4","#41AE76", "#238B45",
    brewer.pal(9, "Greys")[4:8]
  ),order = rev(c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC","Endo","Acinar","Pericyte","EMT_like","Other"))
  )+theme_bw()+
    theme(plot.background = element_blank(),panel.background = element_blank(),panel.grid = element_blank())
  
  Idents(TotalTissue.harmony) = TotalTissue.harmony$David_Tuevson_Clusters
  p2=DimPlot(TotalTissue.harmony,label = TRUE,raster=T,cols = c(
    "#F768A1","#9EBCDA",#"#8C96C6", "#8C6BB1", 
    "#FED976","#FC4E2A", #"#FEB24C", #"#FD8D3C" ,
    "#CCECE6", "#99D8C9","#41AE76", "#238B45",# "#66C2A4",
    brewer.pal(9, "Greys")[4:7]
  ))+theme_bw()+
    theme(plot.background = element_blank(),panel.background = element_blank(),panel.grid = element_blank())
  
  plot = ggpubr::ggarrange(p1,p2,ncol=2)
  pdf(file=paste0("write/Figure_S13a_scRNA-seq_UMAP_",Sys.Date(),".pdf"),width=14,height = 5)
  print(plot)
  dev.off()
  
  meta_df = TotalTissue.harmony@meta.data
  write.csv(meta_df, "write/Figure_S13a_scRNA_seq_meta.csv")
}

# Figure_S13b, Pearson Corr. Coef. using shared cell types
# upper panel: RNA vs RNA, bottom panel: Protein vs Protein
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2024-03-27.rda")
  obj_list = readRDS(file = "ct_dataset_01.rds")
  
  # scRNA-seq
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
    
  av.exp <- AverageExpression(TotalTissue.harmony.filter)$RNA
  cor.exp <- as.data.frame(cor(av.exp))

  # RNA vs RNA
  p1=ggcorr(cor_matrix =cor(av.exp)[],data=NULL, midpoint=0.5,label = FALSE,
         label_round=2,label_size = 4,size=5,low="steelblue",mid="white",high="red",
         legend.position = "bottom",limits=c(0,1),label_alpha = TRUE)#+theme(legend.position = "none")#+coord_flip()
  write.csv(cor(av.exp)[],file="write/Figure_S13b_rna_cor.csv")
  pdf(file=paste0("write/Figure_S13b_corr_rna_",Sys.Date(),".pdf"))
  print(p1)
  dev.off()
  
  # protein vs protein
  meta.df = obj_list$meta
  lfq_mean.df = data.frame(row.names = row.names(obj_list$filter))
  for(i in levels(meta.df$CellType)[-2]){
    temp_index = which(meta.df$CellType==i)  
    lfq_mean.df[,i] = rowMeans(obj_list$filter[,temp_index],
                                            na.rm = TRUE)
  }#darkgrey","white","darkgoldenrod1","brown3"
  
  p2=ggcorr(cor_matrix =cor(lfq_mean.df)[],data=NULL, midpoint =0.5,label = FALSE,
         label_round=2,label_size = 4,size=5,low="steelblue",mid="white",high="red",
         legend.position = "bottom",limits=c(0,1),label_alpha = TRUE)#+theme(legend.position = "none")+coord_flip()
  write.csv(cor(lfq_mean.df)[],file="write/Figure_S13b_protein_cor.csv")
  pdf(file=paste0("write/Figure_S13b_corr_protein_",Sys.Date(),".pdf"))
  print(p2)
  dev.off()

}

# cor, Protein vs RNA
{
    rm(list=ls())    
    obj_list = readRDS(file = "ct_dataset_01.rds")
    ct13_order = levels(meta.df$CellType)[-2] # remove pan-CAF
    copy_num.df = obj_list$identified_copynum # using Perseus based on intensity file

    copy_num.df = na.omit(copy_num.df[,obj_list$meta$SampleID])
    dim(copy_num.df)
    
    meta.df = obj_list$meta
    copynum_sum.df = data.frame(row.names = row.names(copy_num.df))
    for(i in ct13_order){
        temp_index = which(meta.df$CellType==i)  
        copynum_sum.df[,i] = rowMeans(copy_num.df[,temp_index],
                                na.rm = TRUE)
    }
            
    TotalTissue_filter.seurat = readRDS(file = "GSE129455_All_Viable_with_subset_filter_2024-03-27.rda")
    Idents(TotalTissue_filter.seurat) = TotalTissue_filter.seurat$Main_subset_labels_V2
    
    output_cor = c()

    for(i in ct13_order){
        temp_RNA.df = subset(TotalTissue_filter.seurat,idents=c(i))
        temp_RNA.df =exp(temp_RNA.df[["RNA"]]@data)-1
        count_sum.df = apply(temp_RNA.df,1,mean) %>% as.data.frame()
        names(count_sum.df) = "RNA"
        count_sum.df$genename = row.names(count_sum.df)
        #merge
        temp_protein.df = copynum_sum.df %>% dplyr::select(i)
        names(temp_protein.df) = "Protein"
        temp_protein.df$pid = row.names(copynum_sum.df)
        anno.df = obj_list$annotation %>% filter(genename!="")
        comparison.df = anno.df %>% merge(temp_protein.df,by.x="pid",by.y="pid",all.x=F,all.y=F)
        comparison.df = comparison.df %>% merge(count_sum.df,by="genename",all.x=F,all.y=F)
        # comparison.df = comparison.df %>% filter("PCC"!=0 )
        comparison.df = comparison.df[-which(comparison.df$Protein==0 |comparison.df$RNA==0),]
        comparison.df$log2_Protein = log2(comparison.df$Protein)
        comparison.df$log2_RNA = log2(comparison.df$RNA)
        
        cor = cor.test(comparison.df$log2_Protein,comparison.df$log2_RNA,method = "pearson",exact=TRUE)
        output_cor = c(output_cor,cor$estimate)

    }
    cor_df = data.frame(celltype = ct13_order, cor = output_cor)
    write.csv(cor_df, "write/Figure_S13b_RNAvsProtein_cor.csv")
}

# Figure_S13c, Heatmap, known markers, Protein vs RNA
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2024-03-27.rda")
  obj_list = readRDS(file = "ct_dataset_01.rds")
  
  # RNA
  scRNA_meta.df = TotalTissue.harmony@meta.data
  scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
  scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
  scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
  scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
  scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
  
  scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
  scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(obj_list$meta$CellType))
  
  scRNA_meta.df$Main_subset_labels_V2=factor(scRNA_meta.df$Main_subset_labels_V2,levels=(
    c("PCC","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC"
    )))

  TotalTissue.harmony@meta.data = scRNA_meta.df
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels_V2
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
  # average expression of RNA
  RNA_average.df = AverageExpression(TotalTissue.harmony,features = deg_genes)$RNA %>% as.data.frame()
  
  # average intensity of Protein
  lfq.df = obj_list$filter
  
  meta.df = obj_list$meta
  lfq_mean.df = data.frame(row.names = row.names(lfq.df))
  for(i in levels(meta.df$CellType)[-2]){
    temp_index = which(meta.df$CellType==i)  
    lfq_mean.df[,i] = rowMeans(lfq.df[,temp_index],
                               na.rm = TRUE)
  }
  
  anno.df = obj_list$annotation[obj_list$annotation$genename %in% deg_genes2,]
  anno.df = anno.df[c(1:39,41,42,44),]
  lfq_mean.df = lfq_mean.df[anno.df$pid,]
  row.names(lfq_mean.df) = anno.df$genename
  lfq_mean.df = lfq_mean.df[deg_genes2,]
  
  #combine RNA and Protein
  RNA_average.df$order = seq(from=1,to=83,by=2)
  
  lfq_mean.df$order = seq(from=2,to=84,by=2)
  RNA_average.df = RNA_average.df %>% dplyr::select(names(lfq_mean.df))
  
  plot.df = rbind(RNA_average.df,lfq_mean.df)
  plot.df = plot.df[order(plot.df$order),]
  plot.df$order = NULL
  
  plot_anno.df = data.frame(group=factor(rep(c("RNA","Protein"),42),levels = c("RNA","Protein")))
  row.names(plot_anno.df) = row.names(plot.df)
  
  # begin plot
  p1 = pheatmap::pheatmap(plot.df,cellwidth = 18,cellheight = 6,
                     scale = "row",border_color = "white",
                     cluster_rows=FALSE,cluster_cols = FALSE,angle_col=45,
                     annotation_names_row = F,annotation_names_col = F,legend = T,
                     fontsize = 5,display_numbers = FALSE,
                     #silent = FALSE,
                     annotation_row  = plot_anno.df,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                     # color = colorRampPalette(c("blue","white","red"))(100), 
                     width = 9,height = 5)
  # save plot
  pdf(file=paste0("write/Figure_S13c_RNA-Protein-DEG_Pheatmap_",Sys.Date(),".pdf"),width=5,height = 8)
  print(p1)
  dev.off()

  # save output
  mat = as.matrix(plot.df)
  mat = apply(mat,1,scale) %>% t() %>% as.data.frame()
  names(mat) = names(plot.df)
  write.csv(file="write/Figure_S13c_scRNA-seq_pheatmap.csv",x = mat)
  
}

# comparison of cell-percentage from public scRNA-seq dataset
# NOT used
{
    rm(list=ls())
    TotalTissue.harmony = readRDS(file = "GSE129455_All_Viable_with_subset_2024-03-27.rda")
    obj_list = readRDS(file = "ct_dataset_01.rds")
    
    # pre-process of scRNA-seq with shared cell types
    scRNA_meta.df = TotalTissue.harmony@meta.data
    scRNA_meta.df$Main_subset_labels_V2 = as.character(scRNA_meta.df$Main_subset_labels)
    scRNA_meta.df[scRNA_meta.df$Main_subset_labels_V2 %in% c("Endo","Acinar","Other","Pericyte"),"Main_subset_labels_V2"]="Other"
    scRNA_meta_filter.df = scRNA_meta.df %>% filter(Main_subset_labels_V2 != "Other")
    scRNA_meta_filter.df = scRNA_meta_filter.df %>% select(sample,Main_subset_labels_V2)
    scRNA_meta_filter_CAF.df = scRNA_meta_filter.df %>% filter(Main_subset_labels_V2 %in% c("iCAF","myCAF","apCAF"))
    scRNA_meta_filter_CAF.df$Main_subset_labels_V2 = "CAF"
    
    scRNA_meta_filter.df = rbind(scRNA_meta_filter.df,scRNA_meta_filter_CAF.df)
    scRNA_meta_filter.df$Main_subset_labels_V2 = factor(scRNA_meta_filter.df$Main_subset_labels_V2,levels(obj_list$meta$CellType))
    
    scRNA_number.df = scRNA_meta_filter.df %>% group_by(sample,Main_subset_labels_V2) %>% summarize(value=n())
    scRNA_number.df$omic = "scRNA"
    names(scRNA_number.df) = c("Sample","CellType","Value","Omic")
    
    plot.data =scRNA_number.df
    plot.data$CellType = factor(plot.data$CellType,levels=levels(obj_list$meta$CellType))

    # plot the cell-perc with 4 kpc mice
    p1 = ggplot(plot.data, aes(fill=CellType, y=Value, x=Sample)) + 
        geom_bar(position="fill", stat="identity") + 
        theme_classic() + facet_wrap(.~Omic) + xlab("")+ylab("") +
        scale_fill_manual(values = levels(obj_list$meta$CellType_Color)) +
        theme(plot.background = element_blank(),panel.background = element_blank(),strip.background = element_blank(),
            legend.background = element_blank(),text = element_text(size=12),panel.grid = element_blank())#+
    pdf(file=paste0("write/A_CellPerc_RNA_",Sys.Date(),".pdf"),width=7,height = 5)
    print(p1)
    dev.off()
    
    write.csv(plot.data,file=paste0("write/A_CellPerc_RNA_",Sys.Date(),".csv"))
  

}
