library(readxl)
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

setwd("")
# Figure_S12c, heat map of the Pearson Corr. Coef. of all samples
{
  rm(list=ls())
  Obj.list = readRDS(file = "ct_dataset_01.rds")
  
  temp.df = Obj.list$Mini_impute_log2
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
  pdf(file = paste0("Figure_S12c_Corr_Pheatmap_", Sys.Date(),".pdf"), width=7,heigh=6)
  print(p1)
  dev.off()
  
  write.csv(temp.cor, file=paste0("Figure_S12c_Corr_Pheatmap_", Sys.Date(),".csv"))
  
}

# Figure_S12d, heat map for visualization of known markers
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
  
  known_markers_lfq.df = Obj.list$Mini_impute_log2[known_markers.df$pid,]## all samples,markers
  
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
  
  pdf(file = paste0("Figure_12d_Selected_Known_Markers_Pheatmap_",Sys.Date(),".pdf"),
      width=7,heigh=6)
  print(p1)
  dev.off()
  
  write.csv(known_marker_lfq_mean.df, file=paste0("Fig_S12d_Selected_Known_Markers_Pheatmap_",Sys.Date(),".csv",sep=""),row.names = T)
  
}
