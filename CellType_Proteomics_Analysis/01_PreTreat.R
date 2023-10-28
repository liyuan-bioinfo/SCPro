'
---
title: "PreTreat"
author: "Li-Yuan"
date: "2022-10-10"
output: github_document
---
'
## 数据前处理
## 基于面向对象的思想，将FACS-Sorting搜库的数据处理成不同的对象，处理的数据包含**pid**,**intensity.xx**, **LFQ.intensity.xx**等列，前处理的步骤包含**filter**,**normalization**,**imputation**。

## Loading packages and setting global params
rm(list=ls())
library(dplyr)
library(RColorBrewer)

setwd("E:\\XYF（许燕芬）\\FACS_SISPROT\\workspace\\20221128_final/")
options(stringsAsFactors=F)
options(encoding = "UTF-8")



## Loading raw data and experiment data
{
  rm(list=ls())
  ## create list object
  Obj.list = list()
  Obj.list$raw = read.csv(file="../resource/20220612_kpc_cells_msfragger_combined_protein(1).tsv",
                                      fill = FALSE,sep="\t",check.names = F)
  # head(names(Obj.list$raw))
  row.names(Obj.list$raw) = Obj.list$raw$`Protein ID`
  
  Obj.list$raw = Obj.list$raw %>% filter(Organism=="Mus musculus OX=10090") %>% filter(`Indistinguishable Proteins`=="")#7946
  
  ## load experiment data 
  temp_meta.df = read.csv(file="../resource/experiment_design.csv",sep=",",fill=F,check.names = F,header = T)
  
  ## remove NK for following analysis
  temp_meta.df = temp_meta.df[-which(temp_meta.df$CellType == "NK"),]
  temp_meta.df$CellType = factor(temp_meta.df$CellType,
                                 level=c("PCC","CAF","iCAF","myCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC"))
  temp_meta.df$CellLineage = factor(temp_meta.df$CellLineage,
                              levels=c("Cancer Cell", "Fibroblast","Lymphoid Cell", "Myeloid Cell"))
  
  temp_meta.df = temp_meta.df[order(temp_meta.df$CellType),]
  ## set color of different group
  # display.brewer.all(type = "seq")
  # display.brewer.pal(9, "Greys")
  # display.brewer.pal(9, "BuGn")

  # brewer.pal(9, "RdPu")[5] #PCC
  # brewer.pal(9, "YlGnBu")[3:6] #CAF
  # brewer.pal(9, "YlOrRd")[3:6] #Lym
  # brewer.pal(9, "BuGn")[3:7] #Mye
  temp_meta_color.df = temp_meta.df %>% select(CellType) %>% unique()
  temp_meta_color.df$CellType_Color = c(  
    brewer.pal(9, "RdPu")[5], #PCC
    brewer.pal(9, "BuPu")[3:6],#CAF
    brewer.pal(9, "YlOrRd")[3:6],#Lymphoid
    brewer.pal(9, "BuGn")[3:7] )
  temp_meta_color.df$CellType_Color = factor(temp_meta_color.df$CellType_Color,
                                   levels=temp_meta_color.df$CellType_Color)
  
   temp_meta_color.df$CellLineage_Color = c(
    brewer.pal(9, "RdPu")[5], #PCC
    rep(brewer.pal(9, "BuPu")[7],4),#CAF
    rep(brewer.pal(9, "YlOrRd")[7],4),#Lymphoid
    rep(brewer.pal(9, "BuGn")[8],5)
  )
  temp_meta_color.df$CellLineage_Color = factor(temp_meta_color.df$CellLineage_Color, levels=unique(temp_meta_color.df$CellLineage_Color))
  
  temp_meta.df = temp_meta.df %>% merge(temp_meta_color.df,by="CellType")
  temp_meta.df = temp_meta.df[order(temp_meta.df$CellType),]
  Obj.list$meta = temp_meta.df
  
  Obj.list$annotation = Obj.list$raw %>% select("Protein ID","Gene")
  names(Obj.list$annotation) = c("pid","genename")
  
  ## load surface protein as dataframe
  temp.list = readRDS("../resource/MMU.pm.protein.database_4518_20220626.rda")
  Obj.list$uniprot_mouse_loci = temp.list$combined_pm.db
}

## Pre-treat
{
  ## identified pg
  temp_meta.df = read.csv(file="../resource/intensity_experiment_design.csv",sep=",",fill=F,check.names = F,header = T)
  temp_meta.df = temp_meta.df[-which(temp_meta.df$CellType == "NK"),]  
  temp_meta.df$CellType = factor(temp_meta.df$CellType,
                                 level=c("PCC","CAF","iCAF","myCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC"))
  
  temp_meta.df = temp_meta.df[order(temp_meta.df$CellType),]
  temp.df = Obj.list$raw[,temp_meta.df$Raw_SampleID] 
  temp.df = temp.df[apply(temp.df,1,sum)!=0,]
  names(temp.df) = temp_meta.df$SampleID
  Obj.list$identified = temp.df
  
  ## quantified pg
  temp.df = Obj.list$raw[,Obj.list$meta$Raw_SampleID] 
  temp.df = temp.df[apply(temp.df,1,sum)!=0,]
  names(temp.df) = temp_meta.df$SampleID
  Obj.list$quantified = temp.df
}

## filter
{
  #> filter, quantified pg
  temp.valid = c()
  for(i in 1:dim(Obj.list$quantified)[1]){
    temp.df = Obj.list$quantified[i,]
    for(j in levels(temp_meta.df$CellType)){
      temp.id = temp_meta.df[temp_meta.df$CellType %in% j,"SampleID"]
      if(length(which(temp.df[,temp.id]>0))>=2){
        temp.valid = c(temp.valid,i)
        break
      }
    }
  }
  Obj.list$filter = Obj.list$quantified[temp.valid,] #5900
  
}

## Impute
{
## impute 1: no imputate
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
}

## save RDS for following analysis
{
  setwd("../20221128_final/01-pre_treat/")
  temp.df = Obj.list$raw
  temp.df$genename = Obj.list$annotation[row.names(temp.df),"genename"]
  write.csv(x = temp.df,
            file=paste("SupplementTable/Raw_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
    
  temp.df = Obj.list$filter
  temp.df$genename = Obj.list$annotation[row.names(temp.df),"genename"]
  write.csv(x = temp.df,
            file=paste("SupplementTable/Quantified_Filter_LFQ_Intentisty_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
  temp.df = Obj.list$quantified
  temp.df$genename = Obj.list$annotation[row.names(temp.df),"genename"]
  write.csv(x = temp.df,
            file=paste("SupplementTable/Quantified_LFQ_Intentisty_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
  
  temp.df = Obj.list$identified
  temp.df$genename = Obj.list$annotation[row.names(temp.df),"genename"]
  write.csv(x = temp.df,
            file=paste("SupplementTable/Identified_Intentisty_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
  temp.df = Obj.list$no_imputate_log2
  temp.df$genename = Obj.list$annotation[row.names(temp.df),"genename"]
  write.csv(x = temp.df,
            file=paste("SupplementTable/Quantified_log2_LFQ_Intentisty_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
  temp.df = Obj.list$mann_imputate_log2
  temp.df$genename = Obj.list$annotation[row.names(temp.df),"genename"]
  write.csv(x = temp.df,
            file=paste("SupplementTable/Quantified_Imputated_Mann_LFQ_Intentisty_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
  temp.df = Obj.list$uniprot_mouse_loci
  write.csv(x = temp.df,
            file=paste("SupplementTable/Uniport_mouse_loci_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""))
  
  ##> for copynum
  temp.df = Obj.list$identified
  mmu_mass.df = read.table(file="../resource/uniprot-download_true_fields_accession_2Cmass_format_tsv_query__28mo-2022.10.10-07.07.18.89.tsv",header = T)
  temp.df$pid = row.names(temp.df)
  temp_anno.df = temp.df %>% merge(mmu_mass.df,by.x="pid",by.y="Entry",all.x=F,all.y=F)
  write.csv(x = temp_anno.df,
            file=paste("SupplementTable/Identified_Intentisty_Mass_",dim(temp_anno.df)[1],"_",Sys.Date(),".csv",sep=""),quote = F,row.names = F)
  ## after perseus
  {
    temp.df = read.table(file="../resource/intensity_copynum_20221010.txt",check.names = T,fill = T,header = T,sep="\t")
    row.names(temp.df) = temp.df$pid
    temp.df = temp.df[,72:140]
    names(temp.df) = gsub(names(temp.df),pattern = "^Copy.number.",replacement = "")
    write.csv(x = temp.df,
              file=paste("01-pre_treat/SupplementTable/Identified_Intentisty_CopyNum_",dim(temp.df)[1],"_",Sys.Date(),".csv",sep=""),quote = F,row.names = T)
    Obj.list$identified_copynum = temp.df
  }
  
  
  temp_meta.df = Obj.list$meta
  write.csv(x = temp_meta.df,
            file=paste("SupplementTable/Meta_",dim(temp_meta.df)[1],"_",Sys.Date(),".csv",sep=""))
  saveRDS(Obj.list,file = "../FACS.MSFragger.Obj.20221128.rds")
  
}

## membrane annotation V2
library(DOSE)
library(org.Mm.eg.db)
library(clusterProfiler)

## Protein Location Annotation
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  mem.df = Obj.list$annotation[Obj.list$annotation$pid %in%Obj.list$uniprot_mouse_loci$Entry,] %>% unique()
  bg.genename = Obj.list$annotation[Obj.list$annotation$pid %in% row.names(Obj.list$identified),"genename"] %>% unique()
  ##> pid to gene id
  
  genename.id = bitr(unique(mem.df$genename), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  bg.id = bitr(bg.genename, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
  formula_res <- enrichGO(gene = genename.id$ENTREZID, keyType="ENTREZID",
                          OrgDb = org.Mm.eg.db,
                          ont = "CC", pAdjustMethod = "BH",
                          pvalueCutoff = 1,qvalueCutoff = 1,
                          minGSSize = 10,maxGSSize = 5000,
                          universe = bg.id$ENTREZID,
                          readable=TRUE
  )
  
  ##> save data
  write.csv(formula_res@result,file=paste0("01-pre_treat\\SupplementTable\\Identified_Membrane_Specific_GO-CC_geneid_",Sys.Date(),".csv"))
  
  # pdf(file=paste0("Bar_Identified_Mem_GO-CC_genename_",Sys.Date(),".pdf"),
  #     width=8,height=4)
  # barplot(formula_res,color="pvalue",showCategory = 5,label_format = 35,font.size = 12)
  # dev.off()
  
  # pdf(file=paste0("GOgraph_Identified_Mem_GO-CC_genename_",Sys.Date(),".pdf"))
  # plotGOgraph(formula_res)
  # dev.off()
  
  ##> filter surface protein 
  pm_items = c("GO:0005886")
  data.df = formula_res@result 
  data.df = data.df[data.df$ID %in% pm_items,c("ID","Description","pvalue","geneID","Count")]
  
  pm.genename =c()
  for(i in data.df$geneID){
    pm.genename = c(pm.genename,strsplit(i,split = "/")[[1]])
    
  }
  pm.genename = unique(pm.genename) #705
  pm.df = Obj.list$annotation[Obj.list$annotation$genename %in% pm.genename,] #712
  pm.df$loci="pm"
  mem.df$loci="membrane"
  
  all_mem.df = rbind(pm.df,mem.df) %>% select("pid","loci") %>% unique()
  
  write.csv(all_mem.df, file=paste0("01-pre_treat\\SupplementTable\\Identified_Mem_Protein_Go-CC_",dim(all_mem.df)[1],"_",
                                    Sys.Date(),".csv"))
  Obj.list$GO_mouse_pm = all_mem.df
  saveRDS(Obj.list,file = "../FACS.MSFragger.Obj.20221128.rds")
  
  
} 

## PM Protein Function Annotation
{
  protein_anno.df = read.delim(file="../../resource/OmnipathR.interactions_human_48570protein_20210712.txt",header = T)
  protein_functional.df = protein_anno.df %>% filter(aspect == "functional")%>% 
    filter(entity_type=="protein") %>% select(genesymbol,category,parent) %>% 
    filter(genesymbol!="") %>% unique() #21312
  
  Obj.list = readRDS("../FACS.MSFragger.Obj.20221010.rds")
  pm.df = Obj.list$GO_mouse_pm %>% filter(loci=="pm") %>% unique()#712
  surface_protein_mmu.df = Obj.list$annotation[pm.df$pid,]
  surface_protein_mmu.df$genename_hsa = toupper(surface_protein_mmu.df$genename)
  
  surface_category.df = surface_protein_mmu.df %>% merge(protein_functional.df,
                                                         by.x="genename_hsa",by.y="genesymbol",all.x=T,all.y=F)
  
  # write.csv(file=paste0("Category_surface_protein_",Sys.Date(),".csv"),x = surface_category.df)
  GPCR.df = surface_category.df %>% filter(category %in% c("gpcr")) %>% mutate(class="GPCR") %>% select(pid,genename,class) %>% unique()
  ion_channel.df = surface_category.df %>% filter(category %in% c("ion_channel")) %>% mutate(class="ion_channel") %>% select(pid,genename,class) %>% unique()
  cell_adhesion.df = surface_category.df %>% filter(category %in% c("cell_adhesion")) %>% mutate(class="cell_adhesion") %>% select(pid,genename,class) %>% unique()
  transporter.df = surface_category.df %>% filter(category %in% c("transporter")) %>% mutate(class="transporter") %>% select(pid,genename,class) %>% unique()
  RTK.df = surface_category.df %>% filter(category %in% c("cytokine_rtk","receptor_tyrosine_kinases_rtk")) %>% mutate(class="RTK") %>% select(pid,genename,class) %>% unique()
  integrin.df = surface_category.df %>% filter(category %in% c("integrin")) %>% mutate(class="integrin") %>% select(pid,genename,class) %>% unique()
  
  # surface_category_filter.df = surface_category.df %>% filter(category %in% c("gpcr","ion_channel","transporter","cell_adhesion","cytokine_rtk","receptor_tyrosine_kinases_rtk","cytokine_receptor","integrin"))
  
  immune_checkpoint.df = read.delim(file="../../resource/database_20220929.txt",sep="\t",quote = "",check.names = F,header = T,fill = TRUE)
  immune_checkpoint.df = surface_category.df[surface_category.df$genename %in%immune_checkpoint.df$`Rec_Genename(mmu)`,] %>% 
    mutate(class="immune_checkpoint") %>% select(pid,genename,class) %>% unique()
    
  pm_function.df = rbind(GPCR.df,ion_channel.df,cell_adhesion.df,transporter.df,RTK.df,integrin.df,immune_checkpoint.df)#293
  out.df = Obj.list$annotation[pm.df$pid,] %>% merge(pm_function.df,by="pid",all.x=T) 
  out.df$genename.y = NULL
  out.df$class[is.na(out.df$class)]="other"
  
  Obj.list$surface_class = out.df
  
  write.csv(out.df, file=paste0("../01-pre_treat/SupplementTable/Identified_PM_class_",dim(out.df)[1],"_",
                                    Sys.Date(),".csv"))
  saveRDS(Obj.list,file = "../FACS.MSFragger.Obj.20221010.rds")
  
}

##--------
#NOT use
## venn plot membrane db
{
  library(VennDiagram)
  rm(list=ls())
  Obj.list = readRDS(file = "FACS.MSFragger.Obj.20220905.rds")
  db.df = Obj.list$uniprot_mouse_loci
  venn.plot=venn.diagram(
    x = list(
      "Uniprot"= db.df %>% filter(!is.na(uniprot_db)) %>% select(Entry) %>% unlist(),
      "DeepTMHMM"=db.df %>% filter(!is.na(Deep_TMHMM)) %>% select(Entry) %>% unlist(),
      "Phobius"= db.df %>% filter(!is.na(phobius.db)) %>% select(Entry) %>% unlist(),
      "Literature"=db.df %>% filter(!is.na(Surface.PNAS.db)) %>% select(Entry) %>% unlist()
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
    fill = c(alpha("#F768A1",0.3), alpha('#88419D',0.3), alpha('#E31A1C',0.3),alpha('#006D2C',0.3)),
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1,
    cat.default.pos = "outer",
    cat.pos = c(0, 50, 100,200),
    cat.dist = c(0.055, 0.055, 0.085,0.085),
    cat.fontfamily = "sans"
    # rotation = 1
  )
  
  pdf(file=paste("Fig5/MembraneDB_VennPlot_",Sys.Date(),".pdf",sep=""),width = 4,height = 5)  
  grid.draw(venn.plot)
  dev.off()
  
}