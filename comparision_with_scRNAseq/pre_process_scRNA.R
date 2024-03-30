#'@desc workflow for pre-precess of scRNA-seq dataset from GSE129455
#'@author Yuan
#===============================================================================
#           Part I, Basic Analysis with scRNA-seq data using Seurat            #
#================================================================================
rm(list=ls())
library(Seurat)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(Orthology.eg.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
# devtools::install_github("immunogenomics/harmony")
library(harmony)    

options(stringsAsFactors = F)
options(encoding = "UTF-8")

setwd("")
{
  rm(list=ls())
  kpc.data <- read.csv("GSE129455_All_Viable_expression.csv",header = TRUE)
  saveRDS(kpc.data,file="GSE129455_All_Viable_expression.rda")  
}

# scRNA-seq pre-treat
#'@input GSE129455_All_Viable_expression.rda
#'@desc rename genename; basic pre-process of scRNA-seq [normalise and scaled]
{
  temp.KPC.raw = readRDS("GSE129455_All_Viable_expression.rda")
  row.names(temp.KPC.raw)=temp.KPC.raw$X
  temp.KPC.raw$X = NULL
  
  KPC.Seurat <- CreateSeuratObject(counts = temp.KPC.raw, min.cells = 1, min.features = 1)
  temp.meta = KPC.Seurat@meta.data
  temp.meta$orig.ident = row.names(temp.meta)
  temp.meta$barcode = gsub(x = temp.meta$orig.ident,pattern = "\\.KPC.$",replacement = "")
  temp.meta$sample = gsub(x = temp.meta$orig.ident,pattern = "^[ATCG]+\\.",replacement = "")
  
  KPC.Seurat@meta.data = temp.meta
  Idents(KPC.Seurat) = KPC.Seurat$orig.ident
  
  ## Ensembl to Official symbol
  #BiocManager::install("org.Mmu.eg.db")

  ids=AnnotationDbi::select(org.Mm.eg.db, keys=rownames(KPC.Seurat),
             columns=c('ENSEMBL', 'SYMBOL'),
             keytype= 'ENSEMBL')
  
  dim(ids)#13885
  ids=na.omit(ids)#13551
  ids = ids[!duplicated(ids$SYMBOL),]
  ids=ids[!duplicated(ids$ENSEMBL),]
  pos=match(ids$ENSEMBL,rownames(KPC.Seurat))
  KPC.Seurat=KPC.Seurat[pos,]
  
  # RenameGenesSeurat
  RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
    print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
    RNA <- obj@assays$RNA
    
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
      if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
    } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
    obj@assays$RNA <- RNA
    return(obj)
  }
  KPC.Seurat.rename = RenameGenesSeurat(obj = KPC.Seurat, newnames = ids$SYMBOL)

  KPC.Seurat.rename[["RNA"]]@meta.features <- data.frame(row.names = rownames(KPC.Seurat.rename[["RNA"]]))
  #find variable gene
  KPC.Seurat.rename<- FindVariableFeatures(object = KPC.Seurat.rename, mean.function = ExpMean, dispersion.function = LogVMR)

  # regression with Cell Cycle Genes from human
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  temp.egs <- mapIds(org.Hs.eg.db,s.genes, "ENTREZID","SYMBOL")
  mapped <- AnnotationDbi::select(Orthology.eg.db, temp.egs, "Mus.musculus","Homo.sapiens")
  mapped$Mus <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
  # mapped = mapped %>% filter(!is.na(Homo.sapiens))
  temp.s.genes.mus = unique(as.character(unlist(mapped$Mus)))
  mapped[is.na(mapped$Homo.sapiens),]
  temp.s.genes.mus = c(temp.s.genes.mus, "Cenpu")
  
  temp.egs <- mapIds(org.Hs.eg.db,g2m.genes, "ENTREZID","SYMBOL")
  mapped <- AnnotationDbi::select(Orthology.eg.db, temp.egs, "Mus.musculus","Homo.sapiens")
  mapped$Mus <- mapIds(org.Mm.eg.db, as.character(mapped$Mus.musculus), "SYMBOL", "ENTREZID")
  temp.g2m.genes.mus = unique(as.character(unlist(mapped$Mus)))
  # mapped[is.na(mapped$Homo.sapiens),]
  temp.g2m.genes.mus = c(temp.g2m.genes.mus, "Pimreg","Jpt1")
  
  #Calculate Cell Cycle Score (S-G2M Difference)
  KPC.Seurat.rename<- CellCycleScoring(KPC.Seurat.rename, s.features = temp.s.genes.mus,
                                       g2m.features = temp.g2m.genes.mus, set.ident = TRUE)
  KPC.Seurat.rename$CC.Difference <- KPC.Seurat.rename$S.Score - KPC.Seurat.rename$G2M.Score
  
  # Scale
  KPC.Seurat.rename.scale<- ScaleData(object = KPC.Seurat.rename,vars.to.regress = c("nCount_RNA","CC.Difference"))  
  saveRDS(KPC.Seurat.rename.scale,file="GSE129455_All_Viable_expression_scale_20220913.rda")

}

#===============================================================================
#                Part II, correct batch effect using Harmony                   #
#================================================================================
# remove batch effect using Harmony
{  
  rm(list=ls())
  KPC.Seurat = readRDS(file="GSE129455_All_Viable_expression_scale_20220913.rda")
  
  KPC.Seurat = Seurat::RunPCA(KPC.Seurat,features=VariableFeatures(object = KPC.Seurat))
  stdev <- KPC.Seurat@reductions$pca@stdev
  var <- stdev^2
  
  EndVar = 0
  
  for(i in 1:length(var)){
    total <- sum(var)
    numerator <- sum(var[1:i])
    expvar <- numerator/total
    if(EndVar == 0){
      if(expvar > 0.9){
        EndVar <- EndVar + 1
        PCNum <- i
      }
    }
  }
  #Confirm #PC's determined explain > 90% of variance
  sum(var[1:PCNum])/ sum(var)
  
  Idents(KPC.Seurat) = KPC.Seurat$orig.ident
  KPC.Seurat = FindNeighbors(KPC.Seurat,dims = 1:PCNum)
  KPC.Seurat = FindClusters(KPC.Seurat,resolution = 0.5)
  
  KPC.Seurat = Seurat::RunUMAP(KPC.Seurat,dims = 1:PCNum)
  DimPlot(object = KPC.Seurat, reduction = "umap", label = FALSE, pt.size = 0.5,group.by = "sample")
  KPC.Seurat[["UMAP_Clusters"]] = Idents(KPC.Seurat)
  
  # DimPlot(object = TotalTissue.scale, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "DiseaseState", ncol=2)
  
  saveRDS(KPC.Seurat,file="GSE129455_All_Viable_expression_before.rda")
  # load(file="test_out/TestTissue.beforeBatchCorr_202204.RData")
  
  # Batch correction with harmony
  Idents(object = KPC.Seurat) <- 'sample'
  options(repr.plot.height = 2.5, repr.plot.width = 6)
  KPC.Seurat.harmony <- KPC.Seurat %>% 
    RunHarmony("sample", plot_convergence = TRUE)
  
  KPC.Seurat.harmony <- FindNeighbors(object = KPC.Seurat.harmony, dims = 1:PCNum, reduction ="harmony")
  KPC.Seurat.harmony <- FindClusters(object = KPC.Seurat.harmony, resolution = 0.5, reduction ="harmony")
  KPC.Seurat.harmony = RunUMAP(KPC.Seurat.harmony,dims=1:PCNum, reduction = "harmony")
  DimPlot(KPC.Seurat.harmony, reduction = "umap",label = FALSE, pt.size=0.5, group.by = "sample")
  KPC.Seurat.harmony[["UMAP_Clusters"]] = Idents(KPC.Seurat.harmony)
  # Idents(object = TotalTissue.harmony) <- 'ID'
  # levels(TotalTissue.harmony)
  
  saveRDS(KPC.Seurat.harmony,file="GSE129455_All_Viable_expression_after.rda")
  
}

#===============================================================================
#     Part III, Annotation of shared cell types with known markers manually     #
#================================================================================
{
  rm(list=ls())
  TotalTissue.harmony = readRDS("GSE129455_All_Viable_expression_after.rda")
  PCC = c("Epcam", "Krt19", "Mmp7", "Krt8", "Krt18","Cdh1")
  Pericyte=c("Rgs5", "Acta2", "Pdgfrb")
  Fibroblast = c("Lum", "Dcn", "Col1a1", "Col3a1", "Fap")
  T_NK =c("Cd4", "Ctla4", "Nkg7","Gzmb","Itga2","Cd3g","Cd8a")
  # NK=c("Nkg7","Gzmb","Itga2","Itga1")
  B_cell = c("Ms4a1", "Cd79a","Cd19","Ly6d")
  
  Macrophage=c("Adgre1","Saa3","C1qc")
  Monocyte=c("Cd14","Treml4","Cx3cr1","Ccr2")
  Neutrophil=c("S100a8","S100a9","G0s2","Ly6g")
  DCs=c("Ccr7","Ccl22")
  
  # Mast=c("Tpsab1", "Cpa3")
  Acinar=c("Ctrb1", "Prss2","Try5")
  Endothelial=c("Cdh5", "Plvap", "Vwf")
  EMT_like = c("S100a6", "Igfbp4","Sparc","Vim","Spp1")
  Myeloid=c("C1qb","Csf1r","Ear2","Lyz1")
  
  ## PCC
  ## 5,7
  PCC = c("Epcam", "Krt19", "Mmp7", "Krt8", "Krt18","Cdh1","Msln")
  DotPlot(TotalTissue.harmony,features = PCC)
  Seurat::FeaturePlot(TotalTissue.harmony,features = PCC)
  VlnPlot(TotalTissue.harmony,features = PCC)
  DimPlot(TotalTissue.harmony,label = TRUE,label.size = 5)
  
  ## Other
  ## 17,16
  marker_list = c(
    "Ctrb1", "Prss2","Try5", # Acinar,17
    "Cdh5", "Plvap", "Vwf" #Endo,16
  )
  
  ## Fib
  marker_list = c(
    "Rgs5", "Acta2", "Pdgfrb", #Pericyte,18
    "Lum", "Dcn", "Col1a1", "Col3a1", "Fap" #Fib,11
  )
  
  ## B
  #15
  marker_list = c(
    "Ms4a1", "Cd79a","Cd19" 
  )
  DotPlot(TotalTissue.harmony,features = marker_list)
  Seurat::FeaturePlot(TotalTissue.harmony,features = marker_list)
  VlnPlot(TotalTissue.harmony,features = marker_list)
  pdf(file=paste0("Dimplot_beforelabel_",Sys.Date(),".pdf"),height = 5,width = 7)
  DimPlot(TotalTissue.harmony,label = TRUE,label.size = 5)
  dev.off()
  
  ## T,12,13
  marker_list = c(
    "Ptprc","Cd3d","Cd3g",
    "Cd4", "Cd8a","Gzmk","Ctla4", "Nkg7","Gzmb","Pdcd1",
    "Il2ra","Foxp3","Tnfrsf18","Klrg1"
  )
  
  ## Mye
  marker_list=c(
    "C1qc","Ear2","Lyz1" #6
  )
  
  ## Neu
  marker_list=c(
    "S100a8","S100a9","G0s2" #3
  )  
  
  ## MAC
  marker_list=c(
    "Cd14","Ly6c1","Ly6c2","Csf1r","Itgam","Fcgr1","Cx3cr1","Treml4", # Mono
    "Adgre1","Saa3","C1qc","Spp1","C1qa","C1qb","Marco" #1,4
    # ,"Ly6g"
  )
  
  ## DC
  #14
  marker_list=c(
    "Fscn1", "Ccr7",
    "Ccl22","Cd274","Ccl5","H2-Ab1","H2-Eb1","Cd74"
  )

  markers = Seurat::FindMarkers(TotalTissue.harmony,ident.1=c(4))  
  write.csv(markers,file="markers_4.csv")
  
  DotPlot(TotalTissue.harmony, features = unique(c(
    "Foxp3","Klrg1","Tnfrsf18")
  )) + RotatedAxis()
  DotPlot(TotalTissue.harmony, features = unique(c(
    "Ptprc",Myeloid,"Itgam","Apoe","Cd68","Spp1","Cd163","Klrg1","Tnfrsf18",
    Macrophage,Monocyte, Neutrophil)#,EMT_like)
  )) + RotatedAxis()
  DotPlot(TotalTissue.harmony, features = unique(c(
    Epithelial,Pericyte,Fibroblast,T_NK,B_cell,
    Macrophage,Monocyte, Neutrophil ,DCs,Acinar,Endothelial)#,EMT_like)
  )) + RotatedAxis()
  df = TotalTissue.harmony@meta.data
  Idents(TotalTissue.harmony) = TotalTissue.harmony$RNA_snn_res.0.5
  current.cluster.ids <- c(1:19)
  current.cluster.ids[c(1,2,3,5,11)]="MAC"
  current.cluster.ids[c(10)]="MO"
  current.cluster.ids[c(9,14,15)]="DC"
  current.cluster.ids[c(13)]="T"
  current.cluster.ids[c(6,8)]="PCC"
  current.cluster.ids[c(4)]="NEU"
  current.cluster.ids[c(16)]="B"
  current.cluster.ids[c(17)]="Endo"
  current.cluster.ids[c(19)]="Pericyte"
  current.cluster.ids[c(18)]="Acinar"
  current.cluster.ids[c(12)]="CAF"
  current.cluster.ids[c(7)]="MYE"
  names(x = current.cluster.ids) <- levels(x = TotalTissue.harmony)
  TotalTissue.harmony <- RenameIdents(object = TotalTissue.harmony, current.cluster.ids)
  TotalTissue.harmony[["Main_cell_labels"]]<- Idents(object = TotalTissue.harmony)
  
  saveRDS(TotalTissue.harmony,file="GSE129455_All_Viable_annotate.rda")
  
  # T cell subset
  {
    T_NK_cell_subset <- subset(TotalTissue.harmony, idents = "T")
    
    T_NK_cell_subset <- FindVariableFeatures(object = T_NK_cell_subset, mean.function = ExpMean, dispersion.function = LogVMR)
    T_NK_cell_subset <- ScaleData(object = T_NK_cell_subset, vars.to.regress = c("nCount_RNA","CC.Difference"))
    T_NK_cell_subset <- RunPCA(object = T_NK_cell_subset, pc.genes = T_NK_cell_subset@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
    
    #Find number of PCs that gives 90% variance
    st_dev <- T_NK_cell_subset@reductions$pca@stdev
    var <- st_dev^2
    sum(var[1:41])/sum(var)
    #Harmony Batch
    #Find neighbors and clusters WITH harmony batch correction
    options(repr.plot.height = 2.5, repr.plot.width = 6)
    T_NK_cell_subset <- T_NK_cell_subset %>% 
      RunHarmony("sample", plot_convergence = TRUE)
    #Find clusters
    T_NK_cell_subset <- FindNeighbors(object = T_NK_cell_subset, dims = 1:41, reduction = "harmony")
    T_NK_cell_subset <- FindClusters(object = T_NK_cell_subset, resolution = 1,  reduction = "harmony")
    #Run UMAP
    T_NK_cell_subset <- RunUMAP(object = T_NK_cell_subset, dims = 1:41, reduction = 'harmony')
    DimPlot(object = T_NK_cell_subset, reduction = "umap", pt.size = 2, label = T)
    
    DotPlot(T_NK_cell_subset, features = unique(c(
      "Cd3g","Cd3e","Cd4","Icos",
      "Cd8a","Cd8b1",
      "Foxp3","Klrg1","Ikzf2","Sigirr","Tnfrsf18","Gzmb","Nkg7","Ctla4")
    )) + RotatedAxis()
    
    
    #Assess Cluster Quality with FindAllMarkers Function
    # T_cell.markers <- FindAllMarkers(T_NK_cell_subset)
    # df = FindAllMarkers(T_NK_cell_subset)
    # write.csv(file = "subset_Tcells.csv",x = df)
    #Label the clusters
    current.cluster.ids <- c(1:5)
    current.cluster.ids[c(1,5)]="Treg"
    current.cluster.ids[c(4)]="T4"
    current.cluster.ids[c(2)]="T8"
    current.cluster.ids[c(3)]="Other"
    names(x = current.cluster.ids) <- levels(x = T_NK_cell_subset)
    T_NK_cell_subset <- RenameIdents(object = T_NK_cell_subset, current.cluster.ids)
    T_NK_cell_subset[["T_NK_cell_subset_labels"]]<- Idents(object = T_NK_cell_subset)
    saveRDS(T_NK_cell_subset,file = paste0("GSE129455_T_NK_Viable_expression_",Sys.Date(),".rda"))

  }
  
  # CAF subset
  {
    CAF_cell_subset <- subset(TotalTissue.harmony, idents = "CAF")
    
    CAF_cell_subset <- FindVariableFeatures(object = CAF_cell_subset, mean.function = ExpMean, dispersion.function = LogVMR)
    CAF_cell_subset <- ScaleData(object = CAF_cell_subset, vars.to.regress = c("nCount_RNA","CC.Difference"))
    CAF_cell_subset <- RunPCA(object = CAF_cell_subset, pc.genes = CAF_cell_subset@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
    
    #Find number of PCs that gives 90% variance
    st_dev <- CAF_cell_subset@reductions$pca@stdev
    var <- st_dev^2
    sum(var[1:41])/sum(var)
    #Harmony Batch
    #Find neighbors and clusters WITH harmony batch correction
    options(repr.plot.height = 2.5, repr.plot.width = 6)
    CAF_cell_subset <- CAF_cell_subset %>% 
      RunHarmony("sample", plot_convergence = TRUE)
    #Find clusters
    CAF_cell_subset <- FindNeighbors(object = CAF_cell_subset, dims = 1:41, reduction = "harmony")
    CAF_cell_subset <- FindClusters(object = CAF_cell_subset, resolution = 0.5,  reduction = "harmony")
    #Run UMAP
    CAF_cell_subset <- RunUMAP(object = CAF_cell_subset, dims = 1:41, reduction = 'harmony')
    DimPlot(object = CAF_cell_subset, reduction = "umap", pt.size = 2, label = T)
    
    DotPlot(CAF_cell_subset, features = unique(c(
      "Dcn","Pdpn",
      "Pdgfrb","Thy1",
      "Col14a1","Ly6c1","Cxcl12",
      "Cd74","Saa3","H2-Ab1")
    )) + RotatedAxis()
    
    
    #Assess Cluster Quality with FindAllMarkers Function
    # T_cell.markers <- FindAllMarkers(T_NK_cell_subset)
    # df = FindAllMarkers(T_NK_cell_subset)
    # write.csv(file = "subset_Tcells.csv",x = df)
    #Label the clusters
    current.cluster.ids <- c(1:5)
    current.cluster.ids[c(1)]="myCAF"
    current.cluster.ids[c(2,3)]="iCAF"
    current.cluster.ids[c(4,5)]="apCAF"
    names(x = current.cluster.ids) <- levels(x = CAF_cell_subset)
    CAF_cell_subset <- RenameIdents(object = CAF_cell_subset, current.cluster.ids)
    CAF_cell_subset[["CAF_cell_subset_labels"]]<- Idents(object = CAF_cell_subset)
    saveRDS(CAF_cell_subset,file = paste0("GSE129455_CAF_Viable_expression_",Sys.Date(),".rda"))
    
  }
  
  # merge subset
  {
    rm(list=ls())  
    TotalTissue.harmony = readRDS(file="GSE129455_All_Viable_annotate.rda")
    T_NK_cell_subset = readRDS(file="GSE129455_T_NK_Viable_expression_2022-09-14.rda")
    CAF_cell_subset = readRDS(file="GSE129455_CAF_Viable_expression_2022-09-14.rda")
    
    old_meta.df = TotalTissue.harmony@meta.data
    old_meta.df$Main_subset_labels = as.character(old_meta.df$Main_cell_labels)
    T_meta.df = T_NK_cell_subset@meta.data %>% dplyr::select(orig.ident,barcode,T_NK_cell_subset_labels)
    CAF_meta.df = CAF_cell_subset@meta.data %>% dplyr::select(orig.ident,barcode,CAF_cell_subset_labels)
    subset_meta.df = data.frame(orig.ident=c(T_meta.df$orig.ident,CAF_meta.df$orig.ident),
                                barcode=c(T_meta.df$barcode,CAF_meta.df$barcode),
                                subset_label=c(T_meta.df$T_NK_cell_subset_labels,CAF_meta.df$CAF_cell_subset_labels))
    old_meta.df[subset_meta.df$orig.ident,"Main_subset_labels"]=as.character(subset_meta.df$subset_label)
    old_meta.df$Main_subset_labels = factor(old_meta.df$Main_subset_labels,
                                            levels=(c("PCC","Endo","Acinar","Pericyte","myCAF","iCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC","Other")))
    TotalTissue.harmony@meta.data = old_meta.df
    
    Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_subset_labels
    DimPlot(TotalTissue.harmony, reduction = "umap",label = TRUE, pt.size=1)
    
    saveRDS(TotalTissue.harmony,file = paste0("GSE129455_All_Viable_with_subset_",Sys.Date(),".rda"))
  }
  
}

# Plot cell-type specific markers of shared cell types
# NOT USED
{
  rm(list=ls())
  TotalTissue.harmony = readRDS(file="GSE129455_All_Viable_annotate.rda")
  Epithelial = c("Epcam", "Krt19", "Krt8", "Krt18")
  Pericyte=c("Rgs5", "Acta2", "Pdgfrb")
  Fibroblast = c("Lum", "Dcn", "Col1a1", "Col3a1")
  T_NK =c("Cd3g","Cd3e","Cd4", "Cd8a","Ctla4","Foxp3", "Nkg7","Gzmb")
  B_cell = c("Ms4a1", "Cd79a","Cd19","Ly6d")
  Macrophage=c("Adgre1","Saa3","C1qc")
  Monocyte=c("Cx3cr1","Ccr2")
  Neutrophil=c("S100a8","S100a9","G0s2")
  DCs=c("Ccr7","Ccl22")
  Myeloid=c("Ear2","Lyz1")
  # Mast=c("Tpsab1", "Cpa3")
  Acinar=c("Ctrb1", "Prss2","Try5")
  Endothelial=c("Cdh5", "Plvap", "Vwf")
  
  Idents(TotalTissue.harmony) = TotalTissue.harmony$Main_cell_labels
  levels(TotalTissue.harmony) = rev(c("Epithelial","Pericyte","Fibroblast",
                                      "T/NK","B","Myeloid","Mac/Mono","Neutrophil",
                                      "DCs","Acinar",
                                      "Endothelial"))
  pdf(file="Dotplot_KPC_known_marker.pdf",height = 6,width = 12)
  p1=DotPlot(TotalTissue.harmony, features=c(Epithelial, Pericyte, Fibroblast,
                                          T_NK, B_cell, Myeloid,
                                          Macrophage,Monocyte,Neutrophil, DCs, Acinar, 
                                          Endothelial))+ RotatedAxis() +
    theme(axis.text = element_text(size = 14),axis.text.x = element_text(size=13)) +
    xlab("Known Markers")+ylab("Cell Types")
  print(p1)
  dev.off()
  
  pdf(file="UMAP_KPC.pdf",height = 8,width = 14)
  p1 = DimPlot(TotalTissue.harmony, reduction = "umap",label = TRUE,label.size = 6, pt.size=1,raster = FALSE)
  print(p1)
  dev.off()
  
  Idents(TotalTissue.harmony) <- factor(Idents(TotalTissue.harmony),
                                          levels= rev(levels(TotalTissue.harmony)))
  pdf(file="Vlnplot_KPC.pdf",height = 15,width = 12)
  p1 = VlnPlot(TotalTissue.harmony,features=c(
    strsplit(paste(
      "Msln,Tspan8,Dhcr24",
      #Epi
      sep=","
    ),split=",")[[1]]
  ),pt.size = 0.000,ncol = 3,combine = TRUE
  ) 
  print(p1)
  dev.off()
}
