#'@time 20240316
#'@author Yuan
#'@desc Bioinformastic analysis for 14 cell type-proteomics. 

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(limma)
library(pheatmap)
setwd("")

# -------------------------------------------------------------------------
#                            I-Basic pre-treat                            #
# -------------------------------------------------------------------------
{    
    # Desc: Data for Perseus
    # Input: ProteinGroups_Quantified file from MS-search software
    {
        rm(list=ls())
        protein_data_name = "FACS_sorting_Raw_7946_2022-11-28.csv"
        meta_data_name = "meta.txt"    
        celltypes = c("PCC","CAF","iCAF","myCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC")

        protein_df = read.delim(protein_data_name,header=T,check.names=FALSE,sep=",") #7946 * 684
        meta_df = read.delim(meta_data_name,header=T) %>% filter(CellType %in% celltypes)# 18 * 3
        
        ## quantified proteins    
        quantified_protein_df = protein_df[,meta_df$Sample] #7946 * 69
        row.names(quantified_protein_df) = protein_df$`Protein ID`
        names(quantified_protein_df) = meta_df$SampleID
        
        quantified_protein_df <- replace(quantified_protein_df, is.na(quantified_protein_df), 0)

        temp_sum = rowSums(quantified_protein_df, na.rm = TRUE)
        quantified_protein_df = quantified_protein_df[which(temp_sum != 0 ), ] #6406 * 69

        ## save Quantified table for Perseus    
        write.table(quantified_protein_df, file = "Quantified_ct14_6406.txt",sep="\t")    

        ## identified proteins
        meta_df$Sample = gsub(meta_df$Sample, pattern = "MaxLFQ ",replacement="")
        identified_protein_df = protein_df[,meta_df$Sample] #7946 * 69
        row.names(identified_protein_df) = protein_df$`Protein ID`
        names(identified_protein_df) = meta_df$SampleID
        
        identified_protein_df <- replace(identified_protein_df, is.na(identified_protein_df), 0)

        temp_sum = rowSums(identified_protein_df, na.rm = TRUE)
        identified_protein_df = identified_protein_df[which(temp_sum != 0 ), ] #7214 * 69
        write.table(identified_protein_df, file = "Identified_ct14_7214.txt",sep="\t")    

        
    }

    # Desc: Create .RDS and save files from Perseus output, including Quantified, Filter and Impute
    # Output: .RDS
    {
        rm(list=ls())
        Obj_list = list() # save RDS
        RDS_file_name = "dataset_04.rds"

        # prepare annotation files
        anno_df = read.delim("FACS_sorting_Raw_7946_2022-11-28.csv",header=T,check.names=FALSE,sep=",") 
        anno_df = anno_df[,c("Protein ID","Gene")] %>% unique()
        names(anno_df) = c("pids","genenames")
        anno_df$pid = gsub(anno_df$pids,pattern=";.*",replacement="")
        anno_df$genename = gsub(anno_df$genenames,pattern=";.*",replacement="")
        row.names(anno_df) = anno_df$pids

        # prepare meta files
        meta_df = read.delim("meta.txt",header=T,check.names=FALSE) 
        head(meta_df)

        # ct6 order
        ct14_order = c("PCC","CAF","iCAF","myCAF","apCAF","T4","T8","Treg","B","MYE","NEU","MO","MAC","DC")
        ct4_order = c("Cancer Cell", "Fibroblast", "Lymphoid Cell", "Myeloid Cell")

        ct14_color = c(brewer.pal(9, "RdPu")[5], #PCC
                        brewer.pal(9, "BuPu")[3:6],#CAF
                        brewer.pal(9, "YlOrRd")[3:6],#Lymphoid
                        brewer.pal(9, "BuGn")[3:7] )
        ct4_color = c(brewer.pal(9, "RdPu")[5], #PCC
                    rep(brewer.pal(9, "BuPu")[7],4),#CAF
                    rep(brewer.pal(9, "YlOrRd")[7],4),#Lymphoid
                    rep(brewer.pal(9, "BuGn")[8],5))

        # prepare Identified file
        identified_df = read.delim("Identified_ct14_7214.txt",header=T,check.names=FALSE) #6406 * 69
        identified_df = identified_df[,meta_df$SampleID]#re-order
        names(identified_df) = meta_df$SampleID #re-name

        # prepare Quantified file
        quantified_df = read.delim("Quantified_ct14_6406.txt",header=T,check.names=FALSE) #6406 * 69
        quantified_df = quantified_df[,meta_df$SampleID]#re-order
        names(quantified_df) = meta_df$SampleID #re-name

        # prepare Filter file
        filter_df = read.delim("Filter_ct14_5900.txt",header=T,check.names=FALSE) #5900 * 69
        row.names(filter_df) = filter_df$Protein
        filter_df = filter_df[,meta_df$SampleID]
        names(filter_df) = meta_df$SampleID

        # prepare Impute file
        impute_df = read.delim("Impute_ct14_5900.txt",header=T,check.names=FALSE) #5900 * 69
        row.names(impute_df) = impute_df$Protein
        impute_df = impute_df[,meta_df$SampleID]
        names(impute_df) = meta_df$SampleID

        # prepare Impute_Median file
        impute_median_df = read.delim("Impute_median_ct14_5900.txt",header=T,check.names=FALSE) #5900 * 14
        row.names(impute_median_df) = impute_median_df$Protein
        impute_median_df = impute_median_df[,ct14_order]

        # set color of ct14 and ct4
        Obj_list$identified_df = identified_df
        Obj_list$quantified_df = quantified_df
        Obj_list$filter_df = filter_df
        Obj_list$impute_df = impute_df
        Obj_list$impute_median_df = impute_median_df
        Obj_list$meta_df = meta_df
        Obj_list$anno_df = anno_df
        Obj_list$ct14_order = ct14_order
        Obj_list$ct4_order = ct4_order
        Obj_list$ct14_color = ct14_color
        Obj_list$ct4_color = ct4_color
        
        saveRDS(Obj_list, file=RDS_file_name)

        # add impute_mean_df
        Obj_list = readRDS(file="dataset_04.rds") 
        meta_df = Obj_list$meta_df
        data_df = Obj_list$impute_df[,meta_df$SampleID]

        impute_mean_df = data.frame(row.names=row.names(data_df))
        for(i in Obj_list$ct14_order){
            ind = meta_df[which(meta_df$CellType == i),"SampleID"]
            temp_ct_df = apply(data_df[ind], 1, mean) %>% as.data.frame()
            impute_mean_df = cbind(impute_mean_df, temp_ct_df)
        }
        names(impute_mean_df) = Obj_list$ct14_order
        # head(impute_mean_df)
        Obj_list$impute_mean_df = impute_mean_df
        saveRDS(Obj_list, file="dataset_04.rds")

        # add copynum_df
        Obj_list = readRDS(file="dataset_04.rds")         
        meta_df = Obj_list$meta_df
        copynum_df = read.delim(file="copynum_ct14_7241.txt",sep="\t",check.names=FALSE)
        row.names(copynum_df) = copynum_df$pid
        copynum_df = copynum_df[,72:140]
        names(copynum_df) = gsub(names(copynum_df), pattern="^Copy number ", replacement="")
        copynum_df = copynum_df[,Obj_list$meta_df$SampleID]

        # add copynum_mean_df
        copynum_mean_df = data.frame(row.names=row.names(copynum_df))
        for(i in Obj_list$ct14_order){
            ind = meta_df[which(meta_df$CellType == i),"SampleID"]
            temp_ct_df = apply(copynum_df[ind], 1, mean) %>% as.data.frame()
            copynum_mean_df = cbind(copynum_mean_df, temp_ct_df)
        }
        names(copynum_mean_df) = Obj_list$ct14_order        
        Obj_list$copynum_mean_df = copynum_mean_df
        Obj_list$copynum_df = copynum_df
        saveRDS(Obj_list, file="dataset_04.rds")

        #for tangram
        names(Obj_list)
        write.csv(Obj_list$meta_df,file="/aaa/zihanwu/yyyli2/projectx_xu_NC2024/raw/Tangram_deconvolution/ct14_filter_df_5900_meta.csv")
        write.csv(Obj_list$filter_df[Obj_list$meta_df$SampleID],file="/aaa/zihanwu/yyyli2/projectx_xu_NC2024/raw/Tangram_deconvolution/ct14_filter_df_5900.csv")


    }

    # Desc: add PM and Mem database
    # Input: uniprot_mouse_loci.csv
    # Output: .RDS
    {
        rm(list=ls())
        # add PM database
        uniprot_mouse_loci = read.csv("uniprot_mouse_loci.csv",header=TRUE)
        Obj_list = readRDS(file="dataset_04.rds") 

        mem_df = Obj_list$anno_df[Obj_list$anno_df$pid %in% uniprot_mouse_loci$Entry,] %>% unique() #1718
        # bg.genename = Obj.list$annotation[Obj.list$annotation$pid %in% row.names(Obj.list$identified),"genename"] %>% unique()
        ##> pid to gene id
        
        genename_df = bitr(unique(mem_df$genename), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
        # bg.id = bitr(bg.genename, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
        formula_res <- enrichGO(gene = genename_df$ENTREZID, keyType="ENTREZID",
                                OrgDb = org.Mm.eg.db,
                                ont = "CC", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,
                                minGSSize = 1,maxGSSize = 10000,
                                # universe = bg.id$ENTREZID,
                                readable=TRUE
        )
        
        ##> save data
        write.csv(formula_res@result,file=paste0("write/Quantified/PM.csv"))
        
        # pdf(file=paste0("Bar_Identified_Mem_GO-CC_genename_",Sys.Date(),".pdf"),
        #     width=8,height=4)
        # barplot(formula_res,color="pvalue",showCategory = 5,label_format = 35,font.size = 12)
        # dev.off()
        
        # pdf(file=paste0("GOgraph_Identified_Mem_GO-CC_genename_",Sys.Date(),".pdf"))
        # plotGOgraph(formula_res)
        # dev.off()
        
        ##> filter surface protein 
        pm_items = c("GO:0005886")
        pm_items_df = formula_res@result %>% filter(ID %in% pm_items) 
        pm_items_df = pm_items_df[pm_items_df$ID %in% pm_items,c("ID","Description","pvalue","geneID","Count")]
        
        pm.genename =c()
        for(i in pm_items_df$geneID){
            pm.genename = c(pm.genename,strsplit(i,split = "/")[[1]])
            
        }
        pm.genename = unique(pm.genename) #779
        pm_df = Obj_list$anno_df[Obj_list$anno_df$genename %in% pm.genename,] #786        
        pm_df$loci="pm" #786
        mem_df$loci="membrane" #1718

        anno_df = Obj_list$anno_df
        anno_df$loci = ""
        anno_df[which(anno_df$pid %in% mem_df$pid),"loci"] = "mem"
        anno_df[which(anno_df$pid %in% pm_df$pid),"loci"] = "pm"

        Obj_list$uniprot_mouse_loci = uniprot_mouse_loci
        Obj_list$anno_df = anno_df
        saveRDS(Obj_list, file="dataset_04.rds")
    }

    # Desc: add category for PM database
    # Input: uniprot_mouse_loci.csv
    # Output: .RDS
    {
        rm(list=ls())
        # add PM database
        uniprot_mouse_loci = read.csv("uniprot_mouse_loci.csv",header=TRUE)
        Obj_list = readRDS(file="dataset_04.rds")

        mem_df = na.omit(Obj_list$anno_df[uniprot_mouse_loci$Entry,]) #1718
        # bg.genename = Obj.list$annotation[Obj.list$annotation$pid %in% row.names(Obj.list$identified),"genename"] %>% unique()
        ##> pid to gene id
        
        genename_df = bitr(unique(mem_df$genename), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
        # bg.id = bitr(bg.genename, fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
        formula_res <- enrichGO(gene = genename_df$ENTREZID, keyType="ENTREZID",
                                OrgDb = org.Mm.eg.db,
                                ont = "MF", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,
                                minGSSize = 1,maxGSSize = 10000,
                                # universe = bg.id$ENTREZID,
                                readable=TRUE
        )
        
        ##> save data
        write.csv(formula_res@result,file=paste0("write/Quantified/PM_MF.csv"))
                
        ## GO-MF anno
        MF_anno_1 = c("cell adhesion molecule binding", # Cell_adhesion
            "integrin binding", # Integrin
            "ion channel activity", # Ion_channel
            "G protein-coupled receptor activity", # GPCR
            "cell adhesion molecule binding", # Cell_adhesion
            "transmembrane receptor protein tyrosine kinase activity", #RTK
            "transporter activity") # Transporter
        GO_MF_df = formula_res@result %>% dplyr::select(c("ID","Description","pvalue","geneID","Count"))
        GO_MF_df = GO_MF_df %>% filter(Description %in% MF_anno_1) #%>% dplyr::select(c("ID","Description","pvalue","geneID","Count"))        
        
        PM_anno_df = data.frame() # save the anno output
        for(i in 1:dim(GO_MF_df)[1]){
            
            
            temp_genes = strsplit(GO_MF_df[i,"geneID"],split = "/")[[1]]

            temp_df = data.frame(genename=temp_genes,desc=GO_MF_df[i,"Description"],id=GO_MF_df[i,"ID"])
            PM_anno_df = rbind(PM_anno_df, temp_df)
            
        }

        immune_checkpoint_df = read.delim("database_20220929.txt",sep="\t",header=TRUE,check.names=FALSE)
        head(immune_checkpoint_df)
        temp_df = merge(genename_df, immune_checkpoint_df,by.x="SYMBOL",by.y="Rec_Genename(mmu)",all=F) #14
        temp_df = temp_df %>% mutate(genename=SYMBOL,desc="Immune_checkpoint",id="manual_curated") %>% dplyr::select(genename,desc,id)
        dim(temp_df)
        PM_anno_df = rbind(PM_anno_df, temp_df)
        dim(PM_anno_df)#631
        table(PM_anno_df$desc)

        Obj_list$uniprot_mouse_PM_MF = PM_anno_df
        saveRDS(Obj_list, file="dataset_04.rds")
    }    

    # Desc: ANOVA analysis of celltype groups, including 14 cell types and 4 cell lineages
    # Input: Obj_list$impute_df; Obj_list$impute_median_df
    # Output: Obj_list$aov_ct14_df
    # for Fig.5e,f;  Fig.6c,d, Extended Data Fig. 10 and Fig. 12
    # NOT USED
    {
        rm(list=ls())
        Obj_list = readRDS(file="dataset_04.rds") 

        meta_df = Obj_list$meta_df        
        meta_df$CellType = factor(meta_df$CellType, levels = Obj_list$ct14_order)

        data_df = Obj_list$impute_df[,meta_df$SampleID]
        data_median_df = Obj_list$impute_median_df[row.names(data_df),Obj_list$ct14_order]#sort

        protein_num = dim(data_df)[1]
        sample_num = dim(data_df)[2]
        # stage-1 ct-enriched proteins
        ## log2 FC, with median        
        dep_df = data.frame()
        for(i in 1:14){
            temp_df = data_median_df[,i] - data_median_df[,1:14]
            temp_pid = names(which(rowSums(temp_df>log2(1)) == 13))
            temp_df2 = data.frame(pids=temp_pid)
            temp_df2$CellType = names(data_median_df)[i]
            temp_df2$log2FC = apply(temp_df[temp_pid,-i],1,median)        
            dep_df = rbind(dep_df,temp_df2)            
        }

        ## pvalue, with each sample
        row.names(dep_df) = dep_df$pids
        dep_df = dep_df[row.names(data_df),]
        Pvalue = c()
        for(i in 1:protein_num){
            pid_df = data_df[i,] %>% t() %>% as.data.frame()
            names(pid_df) = "pid"
            pid_df$SampleId = row.names(pid_df)
            pid_df$CellType = meta_df$CellType
            model = aov(data=pid_df,pid~CellType)
            # temp_pvalue = summary(model)[[1]]$`Pr(>F)`[1]
            sig_df = TukeyHSD(model)$CellType #which ct is enriched
            enriched_ct = dep_df[i,"CellType"]
            enriched_ct_sig_df = sig_df[grep(row.names(sig_df),pattern=paste0("-",enriched_ct,"$|","^",enriched_ct,"-")),] #whethor this ct is sig. with other cts
            temp_pvalue = max(enriched_ct_sig_df[,"p adj"])                                    
            Pvalue = c(Pvalue,temp_pvalue)            
        }        

        dep_df$pvalue = Pvalue
        dep_df$fdr = p.adjust(Pvalue,method = "BH")
        ct14_output_df = merge(dep_df, Obj_list$anno_df, by.x="pids",by.y="pids",all.y=F)
        write.csv(ct14_output_df,"write/Stat/temp_sig_ct14.csv")
        head(output_df)

        # stage-2 cell lineage enriched proteins
        data_df = Obj_list$impute_df[,meta_df$SampleID]#sort
        data_median_df = data.frame(row.names=row.names(data_df))
        for(i in Obj_list$ct4_order){
            temp_index = which(meta_df$CellLineage == i)
            temp_df = apply(data_df[,temp_index],1,median)            
            data_median_df = cbind(data_median_df, temp_df)
        }
        names(data_median_df) = Obj_list$ct4_order

        #log2FC within cell lineage        
        cl_dep_df = data.frame()
        for(i in 1:4){
            temp_df = data_median_df[,i] - data_median_df[,1:4]
            temp_pid = names(which(rowSums(temp_df>log2(1)) == 3))
            temp_df2 = data.frame(pids=temp_pid)
            temp_df2$CellLineage = names(data_median_df)[i]
            temp_df2$log2FC = apply(temp_df[temp_pid,-i],1,median)        
            cl_dep_df = rbind(cl_dep_df,temp_df2)            
        }    
        
        row.names(cl_dep_df) = cl_dep_df$pids
        cl_dep_df = cl_dep_df[row.names(data_df),]#sort

        Pvalue = c()
        for(i in 1:protein_num){
            pid_df = data_df[i,] %>% t() %>% as.data.frame()
            names(pid_df) = "pid"
            pid_df$SampleId = row.names(pid_df)
            pid_df$CellLineage = meta_df$CellLineage
            model = aov(data=pid_df,pid~CellLineage)
            
            # temp_pvalue = summary(model)[[1]]$`Pr(>F)`[1]
            sig_df = TukeyHSD(model)$CellLineage #which ct is enriched
            enriched_ct = cl_dep_df[i,"CellLineage"]
            enriched_ct_sig_df = sig_df[grep(row.names(sig_df),pattern=paste0("-",enriched_ct,"$|","^",enriched_ct,"-")),] #whethor this ct is sig. with other cts

            # temp_pvalue = summary(model)[[1]]$`Pr(>F)`[1]        
            temp_pvalue = max(enriched_ct_sig_df[,"p adj"])                            
            Pvalue = c(Pvalue,temp_pvalue)            
        }        

        cl_dep_df$pvalue = Pvalue
        cl_dep_df$fdr = p.adjust(Pvalue,method = "BH")    

        ct4_output_df = merge(cl_dep_df, Obj_list$anno_df, by.x="pids",by.y="pids",all.y=F)
        write.csv(ct4_output_df,"write/Stat/temp_sig_cl4.csv")

        # output
        Obj_list$aov_ct14_df = ct14_output_df
        Obj_list$aov_ct4_df = ct4_output_df

        saveRDS(Obj_list, file="dataset_04.rds")
    }    
    
    # Desc: LIMMA analysis of celltype groups, including 14 cell types and 4 cell lineages
    # Input: Obj_list$impute_df; Obj_list$impute_median_df
    # Output: Obj_list$limma_ct14_df
    # for Fig.5e,f;  Fig.6c,d, Extended Data Fig. 10 and Fig. 12
    {
        rm(list=ls())        
        Obj_list = readRDS(file="dataset_04.rds") 

        meta_df = Obj_list$meta_df        
        meta_df$CellType = factor(meta_df$CellType, levels = Obj_list$ct14_order)

        data_df = Obj_list$impute_df[,meta_df$SampleID]
        data_median_df = Obj_list$impute_median_df[row.names(data_df),Obj_list$ct14_order]#sort

        protein_num = dim(data_df)[1]
        sample_num = dim(data_df)[2]

        # stage-1 ct-enriched proteins                
        stat_ct14_df = data.frame()
        for(i in 1:length(Obj_list$ct14_order)){
            
            design = model.matrix(~0+CellType, data=meta_df)
            colnames(design) = Obj_list$ct14_order
            ct = Obj_list$ct14_order[i] #this cell type
            contrasts = paste0(ct,"-",Obj_list$ct14_order[-i]) # this ct vs the rest
            contrast.matrix = makeContrasts(contrasts=contrasts, levels=Obj_list$ct14_order)

            # xx <- c("B-A","C-B","C-A")
            # makeContrasts(contrasts=x,levels=c("A","B","C"))
            fit = lmFit(data_df, design)
            fit2 = contrasts.fit(fit, contrast.matrix)
            fit2 = eBayes(fit2)

            pval_df = fit2$p.value
            merged_pvals = apply(pval_df, 1, FUN = function(x) {  
                RecordTest::fisher.method(as.numeric(unlist(x)))$"p.value"[[1]]  
            })            

            log2FC_df = fit2$coefficients
            min_log2FC = apply(log2FC_df, 1, FUN = function(x) {  
                min(as.numeric(unlist(x)))
            })            
            mean_log2FC = apply(log2FC_df, 1, FUN = function(x) {  
                mean(as.numeric(unlist(x)))
            })            
            merged_log2FC = apply(log2FC_df, 1, FUN = function(x) {  
                sum(as.numeric(unlist(x)))
            })            

            temp_merged_df = data.frame(pids=row.names(data_df),sum_pvals=merged_pvals, sum_log2fc=merged_log2FC, min_log2fc=min_log2FC,mean_log2fc=mean_log2FC, CellType=ct)

            stat_ct14_df = rbind(stat_ct14_df, temp_merged_df)
        }        

        stat_ct14_df$fdr = p.adjust(stat_ct14_df$sum_pvals,method = "BH")
        stat_ct14_df = merge(stat_ct14_df, Obj_list$anno_df, by.x="pids",by.y="pids",all.y=F)

        output_stat_ct14_df = stat_ct14_df %>% dplyr::filter(fdr<0.05 & min_log2fc>0)#5623        

        Obj_list$limma_ct14_df = output_stat_ct14_df
        saveRDS(Obj_list, file="dataset_04.rds")

    }

    # Desc: LIMMA analysis of celltype groups, including 14 cell types and 4 cell lineages
    # Input: Obj_list$impute_df; Obj_list$impute_median_df
    # Output: Obj_list$limma_ct4_df
    # for Fig.5e,f;  Fig.6c,d, Extended Data Fig. 10 and Fig. 12
    {
        rm(list=ls())        
        Obj_list = readRDS(file="dataset_04.rds") 
        ct_order = c("PCC","FIB","LYM","MYE")

        meta_df = Obj_list$meta_df      
        meta_df$CellType =  meta_df$CellLineage
        meta_df$CellType = gsub(x=meta_df$CellType,pattern="^Cancer Cell$",replacement="PCC")
        meta_df$CellType = gsub(x=meta_df$CellType,pattern="^Fibroblast",replacement="FIB")
        meta_df$CellType = gsub(x=meta_df$CellType,pattern="^Lymphoid Cell$",replacement="LYM")
        meta_df$CellType = gsub(x=meta_df$CellType,pattern="^Myeloid Cell$",replacement="MYE")
        meta_df$CellType = factor(meta_df$CellType, levels = ct_order)
        meta_df = meta_df %>% arrange(CellType)
        
        data_df = Obj_list$impute_df[,meta_df$SampleID]        

        stat_ct4_df = data.frame()
        for(i in 1:length(ct_order)){
            
            design = model.matrix(~0+CellType, data=meta_df)
            colnames(design) = ct_order
            ct = ct_order[i] #this cell type
            contrasts = paste0(ct,"-",ct_order[-i]) # this ct vs the rest
            contrast.matrix = makeContrasts(contrasts=contrasts, levels=ct_order)
            
            fit = lmFit(data_df, design)
            fit2 = contrasts.fit(fit, contrast.matrix)
            fit2 = eBayes(fit2)

            pval_df = fit2$p.value
            merged_pvals = apply(pval_df, 1, FUN = function(x) {  
                RecordTest::fisher.method(as.numeric(unlist(x)))$"p.value"[[1]]  
            })            

            log2FC_df = fit2$coefficients
            min_log2FC = apply(log2FC_df, 1, FUN = function(x) {  
                min(as.numeric(unlist(x)))
            })            
            mean_log2FC = apply(log2FC_df, 1, FUN = function(x) {  
                mean(as.numeric(unlist(x)))
            })            
            merged_log2FC = apply(log2FC_df, 1, FUN = function(x) {  
                sum(as.numeric(unlist(x)))
            })            

            temp_merged_df = data.frame(pids=row.names(data_df),sum_pvals=merged_pvals, sum_log2fc=merged_log2FC, min_log2fc=min_log2FC,mean_log2fc=mean_log2FC, CellType=ct)

            stat_ct4_df = rbind(stat_ct4_df, temp_merged_df)
        }        

        stat_ct4_df$fdr = p.adjust(stat_ct4_df$sum_pvals,method = "BH")
        stat_ct4_df = merge(stat_ct4_df, Obj_list$anno_df, by.x="pids",by.y="pids",all.y=F)

        output_stat_ct4_df = stat_ct4_df %>% dplyr::filter(fdr<0.05 & min_log2fc>0)#4526        

        write.csv(stat_ct4_df, file="write/Stat/llimma_sig_ct4.csv")

        Obj_list$limma_ct4_df = output_stat_ct4_df
        saveRDS(Obj_list, file="dataset_04.rds")

    }

    # Desc: Identified of novel functional cell-types
    # Input: Obj_list$copynum_mean_df; Obj_list$uniprot_mouse_PM_MF; Obj_list$limma_ct14_df
    # Output: Obj_list$sig_mf_ct14_df
    # for Fig.6a and Fig.6d
    {
        rm(list=ls())        
        Obj_list = readRDS(file="dataset_04.rds") 
        ct_order = c("PCC","FIB","LYM","MYE")

        sig_ct14_df = Obj_list$limma_ct14_df %>% filter(sum_pvals<0.05 & sum_log2fc>1 & min_log2fc>0) #5570
        row.names(sig_ct14_df) = sig_ct14_df$pids

        # anno PM proteins
        pm_df = Obj_list$anno %>% dplyr::filter(loci == "pm")
        sig_ct14_df = na.omit(sig_ct14_df[pm_df$pids,])# 475 PM and Sig. Proteins

        # anno the molecular function for pm proteins
        mf_df = Obj_list$uniprot_mouse_PM_MF
        sig_ct14_df = merge(sig_ct14_df, mf_df, by="genename",all.x=T,all.y=F)
        head(sig_ct14_df)

        # anno copy number
        copynum_mean_df = log2(Obj_list$copynum_mean_df+1)
        copynum_mean_df$pids = row.names(copynum_mean_df)
        copynum_mean_df = copynum_mean_df %>% tidyr::gather(key="CellType",value="copynum",-pids)
                
        log2copynum = c()
        for (i in (1:dim(sig_ct14_df)[1])){
            ct = sig_ct14_df[i,"CellType"]
            pid = sig_ct14_df[i,"pid"]
            temp_cp_df = copynum_mean_df %>% filter(CellType==ct & pids == pid)
            copynum = temp_cp_df$copynum
            if(length(copynum)==0){
                copynum = 0
            }
            
            log2copynum = c(log2copynum, copynum)


        }
        sig_ct14_df$log2copynum = log2copynum
        sig_ct14_df$score = sig_ct14_df$sum_log2fc * sig_ct14_df$log2copynum

        Obj_list$sig_mf_ct14_df = sig_ct14_df

        write.csv(sig_ct14_df, file=paste0("write/Stat/Fig_6a_MF_ct14_",Sys.Date(),".csv"))
        
        saveRDS(Obj_list, file="dataset_04.rds")

    }

}

# -------------------------------------------------------------------------
#                            II-Quality Control                           #
# -------------------------------------------------------------------------
{
    # Desc: Bar plot for quantified protein number
    # Input: Obj_list$quantified_df and identified_df
    # Fig.5b
    {
        rm(list=ls())
        Obj_list = readRDS(file = "dataset_04.rds")
        pm_df = Obj_list$uniprot_mouse_loci
        temp_df = Obj_list$identified_df #%>% t() %>% as.data.frame()
        
        meta_df = Obj_list$meta_df %>% dplyr::select(SampleID,CellType) %>% unique()
        # row.names(plot.data) = plot.data$group
        plot.data = meta_df %>% dplyr::select(CellType) %>% unique()
        row.names(plot.data) = plot.data$CellType
        
        temp_sum = c()
        for(i in Obj_list$ct14_order){
            temp_index = which(meta_df$CellType==i)  
            tt = temp_df[temp_index,]
            ttt = apply(tt,2,function(x){
            return(sum(x,na.rm = TRUE))
            })
            temp_len = length(which(ttt >0))
            temp_sum = c(temp_sum,temp_len)
            
        }
        plot.data$protein_sum =temp_sum
        
        ## quantified
        plot.data2 = meta_df %>% dplyr::select(CellType) %>% unique()
        row.names(plot.data2) = plot.data2$CellType
        temp_df2 = Obj_list$quantified_df %>% t() %>% as.data.frame()
        
        temp_sum = c()
        for(i in Obj_list$ct14_order){
            temp_index = which(meta_df$CellType==i)  
            tt = temp_df2[temp_index,]
            ttt = apply(tt,2,function(x){
            return(sum(x,na.rm = TRUE))
            })
            temp_len = length(which(ttt >0))
            temp_sum = c(temp_sum,temp_len)
            
        }
        plot.data2$protein_sum =temp_sum
        
        #combine
        plot.data["Total","protein_sum"]=dim(Obj_list$identified_df)[1]
        plot.data$CellType = row.names(plot.data)
        plot.data$protein_loca = "Identified"
        
        plot.data2["Total","protein_sum"]=dim(Obj_list$quantified_df)[1]
        plot.data2$CellType = row.names(plot.data2)
        plot.data2$protein_loca = "Quantified"
        
        plot.data3 = plot.data
        plot.data3$protein_sum = plot.data$protein_sum - plot.data2$protein_sum
        plot.data.bar = rbind(plot.data3,plot.data2)
        plot.data.bar$CellType = factor(x = plot.data.bar$CellType,levels = c(Obj_list$ct14_order,"Total"))
        
        temp.p1 <- ggplot(plot.data.bar, aes(x=CellType,y=protein_sum,fill=protein_loca)) + 
            geom_bar(stat = "identity",width=0.8)+
            scale_fill_grey(start = 0.75, end = 0.5) +
            theme_classic()+
            theme(        
            plot.title = element_text(size=16,hjust = 0.5),
            axis.title.y.left  = element_text(size=16),
            axis.text = element_text(size=12),
            axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1,color = c(levels(Obj_list$meta$CellType_Color),"Black")),
            axis.ticks.x = element_blank(),axis.line.x.bottom = element_blank(),
            legend.position=c(0.5,0.8),legend.direction = "horizontal",legend.text = element_text(size=12),
            legend.background = element_blank(),panel.grid = element_blank(),
            plot.background = element_blank(),strip.background = element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank()
            )+
            labs(y="Protein groups(number)",fill="",x="")+
            scale_y_continuous(expand = c(0,0),limits=c(0,8000),breaks=c(0,2000,4000,5000,6000,7000,8000),label=c(0,2000,4000,5000,6000,7000,8000))
            
        
        ## Output  
        plot.data.bar2 = rbind(plot.data,plot.data2)
        write.csv(plot.data.bar2,
                    file = paste0("write/Quantified/Fig5_b_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".csv"),
                    fileEncoding = "UTF-8",row.names = F,quote = TRUE)
        
        pdf(paste0("write/Quantified/Fig5_b_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".pdf"),
                width=6,height = 4)
        print(temp.p1)
        dev.off()
    }
            

    # Desc: Pearson Corr of 14 celltypes
    # Input: Obj_list$quantified_df
    # Output: Pheatmap plot[Fig.5c]; Related table
    {
        rm(list=ls())
        Obj_list = readRDS(file = "dataset_04.rds")        
        meta_df = Obj_list$meta
        # data_df = Obj_list$quantified_df[,meta_df$SampleID]
        data_df = Obj_list$impute_df[,meta_df$SampleID]
        # data_df = log2(data_df+1)
                
        plot_df = cor(data_df[])

        # search the max value
        row_max = c()#0.8422707
        for(i in 1:dim(plot_df)[1]){
            temp_row = plot_df[i,]
            row_max = c(row_max,max(temp_row[which(temp_row!=1)]))            
            
        }

        # replace "1" with max corr
        row_max = max(row_max) # 0.8422707
        for(i in 1:dim(plot_df)[1]){
            temp_row = plot_df[i,]            
            plot_df[i,which(temp_row==1)]=row_max            
        }        

        # set anno_df
        meta_df = meta_df %>% dplyr::select(CellType)
        row.names(meta_df) = Obj_list$meta$SampleID

        # set color of anno_df
        mycolors <- c(Obj_list$ct14_color)
        names(mycolors) <- Obj_list$ct14_order
        mycolors <- list(CellType = mycolors)
        
        p1 = pheatmap(plot_df,cluster_rows = F,cluster_cols = F,width = 7,height = 6,#cellwidth = 3,cellheight = 3,
                        show_rownames = T,show_colnames = F,border_color = NA,
                        annotation_col = meta_df,annotation_names_row = F,annotation_names_col = F,
                        annotation_row = meta_df,annotation_colors = mycolors,
                        annotation_legend = F,fontsize_row = 6,fontsize_col = 6,
                        color = colorRampPalette(c("darkgrey","white","darkgoldenrod1","brown3"))(100),silent = TRUE)


        llim = round(min(plot_df),2) # 0.270
        ulim = round(max(plot_df),2) # 0.842

        # Fig lengends           
        df <- data.frame(matrix(nrow = 10, ncol = 2))
        df[] <- rnorm(20)
        p2=ggplot(data = df, aes(x = X1, y = X2, fill = X2)) + geom_tile() + theme_void()+
            scale_fill_gradientn(colors=colorRampPalette(c("darkgrey","white","darkgoldenrod1","brown3"))(100),
            limits = c(llim, ulim), breaks = round(seq(llim, ulim, length.out=4),2),guide = guide_colorbar(barwidth = 5, barheight = 1, 
                                                                                            ticks = TRUE, ticks.colour="black",ticks.linewidth=1/.pt,direction="horizontal",#"vertical",                                                                                            
                                                                                            label.theme = element_text(size = 8),
                                                                                            draw.ulim = TRUE, 
                                                                                            draw.llim = TRUE, 
                                                                                            draw.separator = TRUE,
                                                                                            
                                                                                            title.position = "top"
                                                                                            )
            )        

        #Plot
        pdf(file = paste0("write/Quantified/Fig5_c_Corr_Pheatmap_Quantified_",
                            Sys.Date(),".pdf"), width=7,heigh=6)
        print(p1)
        print(p2)
        dev.off()
        
        write.csv(cor(data_df[]), file=paste0("write/Quantified/Fig5_c_Corr_Pheatmap_Quantified_",
                                        Sys.Date(),".csv"))
    }


    # Desc: Quantified analysis of 14 celltypes
    # Input: Obj_list$quantified
    # Output: PCA plot[Fig.4d]
    {  
        rm(list=ls())
        Obj_list = readRDS(file = "dataset_04.rds")
        meta_df = Obj_list$meta_df
        
        data_df = Obj_list$quantified_df[,meta_df$SampleID]#select and sort
        data_df = log2(data_df+1) %>% t() %>% as.data.frame()
                
        pca_result <- prcomp(data_df,scale=T,center = T)        
        plot_df = as.data.frame(pca_result$x)

        summ1 <- summary(pca_result)
        xlab1 <- paste0("PC1 (",round(summ1$importance[2,1]*100,2),"%)")
        ylab1 <- paste0("PC2 (",round(summ1$importance[2,2]*100,2),"%)")
                
        plot_df = cbind(plot_df, meta_df)
        plot_df$CellType = factor(plot_df$CellType, levels = Obj_list$ct14_order)        

        p1=ggplot(data = plot_df,aes(x = PC1,y = PC2,color = CellType))+
            stat_ellipse(aes(fill = CellLineage),
                        type = "norm",geom = "polygon",alpha = 0.25,
                        color = NA
            )+ 
            geom_point(size = 3.5)+
            labs(x = xlab1,y = ylab1,color = "",title = "")+
            guides(fill = "none")+
            theme_bw()+
            scale_fill_manual(values = unique(Obj_list$ct4_color))+
            scale_colour_manual(values = unique(Obj_list$ct14_color))+
            theme(plot.background = element_blank(),legend.background = element_blank(),
                panel.background = element_blank(),panel.grid = element_blank(),
                axis.text = element_text(size = 12),axis.title = element_text(size = 16),
                legend.text = element_text(size = 16),
                plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))

        # save plot        
        pdf(file=paste0("write/Quantified/Fig5_d_PCA_CellLineage_Quantified_",Sys.Date(),".pdf"),
            width = 7,height=5)
        print(p1)
        dev.off()
        
        write.csv(x = plot_df,file=paste0("write/Quantified/Fig5_d_PCA_CellLineage_Quantified_",Sys.Date(),".csv"))
    }

    # Desc: Bar plot for identified PM protein number
    # Input: Obj_list$identified_df; Obj_list$uniprot_mouse_loci
    # Fig.5b
    {
        rm(list=ls())
        Obj_list = readRDS(file = "dataset_04.rds")

        meta_df = Obj_list$meta_df %>% dplyr::select(SampleID,CellType) %>% unique()
        pm_anno_df = Obj_list$anno_df %>% dplyr::filter(loci == "pm")
        anno_df = Obj_list$anno_df

        pm_data_df = na.omit(Obj_list$quantified_df[pm_anno_df$pids,meta_df$SampleID])
        data_df = na.omit(Obj_list$quantified_df[anno_df$pids,meta_df$SampleID])
                        
        identified_at_leat_one_ct = c()
        identified_pm_at_leat_one_ct = c()
        for(i in Obj_list$ct14_order){
            temp_index = which(meta_df$CellType==i) # Samples ID
            sum_data_df = apply(data_df[,temp_index],1,function(x){return(sum(x,na.rm = TRUE))}) # sum the cts
            sum_pm_data_df = apply(pm_data_df[,temp_index],1,function(x){return(sum(x,na.rm = TRUE))}) # sum the cts
            identified_at_leat_one_ct = c(identified_at_leat_one_ct,length(which(sum_data_df >0)))            
            identified_pm_at_leat_one_ct = c(identified_pm_at_leat_one_ct,length(which(sum_pm_data_df >0)))            
        }

        count_df = data.frame("identified"=identified_at_leat_one_ct, "identified_pm"=identified_pm_at_leat_one_ct)
        count_df$pm_perc = round(count_df$identified_pm / count_df$identified,4) *100
        # count_df$ct = Obj_list$ct14_order

        # add total
        count_df = count_df %>% t() %>% as.data.frame()
        names(count_df) = Obj_list$ct14_order
        count_df$All = c(dim(data_df)[1], dim(pm_data_df)[1], dim(pm_data_df)[1]/dim(data_df)[1] * 100)
        
        # for plot
        plot_df = count_df %>% t() %>% as.data.frame()
        plot_df$CellType = factor(row.names(plot_df), levels=row.names(plot_df))

        p1 <- ggplot(plot_df, aes(x=CellType)) + 
            geom_bar(aes(y=identified_pm),stat = "identity",width=0.3,fill="grey") +
            geom_line(aes(y=pm_perc*40),stat = "identity",size=1.5,group=1,fill="black") +
            scale_y_continuous(name="identified_pm", sec.axis = sec_axis(~./40,name="pm_perc")) +
            # scale_fill_grey(start = 0.75, end = 0.5) +
            theme_classic()+
            theme(        
                plot.title = element_text(size=16,hjust = 0.5),
                axis.title.y.left  = element_text(size=16),
                axis.text = element_text(size=12),
                axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1,color = c(Obj_list$ct14_color,"Black")),
                axis.ticks.x = element_blank(),axis.line.x.bottom = element_blank(),
                legend.position=c(0.5,0.8),legend.direction = "horizontal",legend.text = element_text(size=12),
                legend.background = element_blank(),panel.grid = element_blank(),
                plot.background = element_blank(),strip.background = element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background = element_blank()
            )+
            labs(y="Protein groups(number)",fill="",x="")#+
            # scale_y_continuous(expand = c(0,0),limits=c(0,600),breaks=c(0,200,400,600),label=c(0,200, 400, 600))
            
        
        pdf(paste0("write/Quantified/Fig6_b_PM_Barplot_",Sys.Date(),".pdf"),
                width=6,height = 4)
        print(p1)
        dev.off()

    }    

    # Desc: Pie plot for molecular function of plasma proteins
}

# -------------------------------------------------------------------------
#                            III-Stat Analysis                            #
# -------------------------------------------------------------------------
{   
    # Desc: Pheatmap of sig.protein panel of 4 cl
    # Input: Obj_list$limma_ct4_df;
    # Output: Pheatmap[Fig.5e]; Related table
    {
        rm(list=ls())
        Obj_list = readRDS(file="dataset_04.rds") 

        meta_df = Obj_list$meta_df        
        meta_df$CellType = factor(meta_df$CellType, levels = Obj_list$ct14_order)
        data_df = Obj_list$impute_df[,meta_df$SampleID]

        ct14_dep_df = Obj_list$limma_ct4_df %>% dplyr::filter(fdr<0.05 & sum_log2fc>1)#1636
        ct_order = c("PCC","FIB","LYM","MYE")
        ct14_dep_df$CellType = factor(ct14_dep_df$CellType, levels = ct_order)
        ct14_dep_df = ct14_dep_df %>% dplyr::arrange(CellType)
        
        plot_df = data_df[ct14_dep_df$pids,]                                        

        ## annotation
        annotation_col = data.frame(CellType = meta_df$CellType)
        row.names(annotation_col) = meta_df$SampleID
        anno_colors = list(CellType = unique(Obj_list$ct14_color))
        names(anno_colors$CellType) = Obj_list$ct14_order
                
        # 使用apply函数对data.frame的每一行应用缩放函数  
        df_scaled_rows <- t(apply(plot_df, 1, scale))  
        
        # 计算99%分位数
        ulim <- quantile(df_scaled_rows, 0.99)
        llim <- quantile(df_scaled_rows, 0.01)

        # 将超过99%分位数的值替换为99%分位数
        df_scaled_rows <- ifelse(df_scaled_rows > ulim, ulim, df_scaled_rows)
        df_scaled_rows <- ifelse(df_scaled_rows < llim, llim, df_scaled_rows)

        # 将结果转换回data.frame格式  
        df_scaled_rows <- as.data.frame(df_scaled_rows)
        names(df_scaled_rows) = names(plot_df)

        length(which(df_scaled_rows>0)) #62
        length(which(df_scaled_rows<0)) #71
        # 创建颜色映射
        my_palette <- c(
            c(colorRampPalette(c("navy", "white"))(50)),
            c(colorRampPalette(c("white", "red"))(50))
        )
        # 将0设置为断点
        my_breaks <- unique(c(seq(min(df_scaled_rows), 0, length.out = 50), seq(0, max(df_scaled_rows), length.out = 50)))

        p1=pheatmap::pheatmap(df_scaled_rows,cluster_rows = F,show_rownames = F,cluster_cols = F,
                scale = "none",border_color = NA,legend = T,cellwidth = 5,
                # color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                color = my_palette,breaks=my_breaks,#colorRampPalette(c("navy", "white", "red"))(100),
                annotation_col = annotation_col,show_colnames = F,legend_labels = F,annotation_legend = F,labels_col = NA,labels_row = NA,
                annotation_names_col = FALSE,annotation_colors = anno_colors,silent = TRUE
        )

        # Fig lengends           
        df <- data.frame(matrix(nrow = 10, ncol = 2))
        df[] <- rnorm(20)
        p2=ggplot(data = df, aes(x = X1, y = X2, fill = X2)) + geom_tile() + theme_void()+
            scale_fill_gradientn(colors=my_palette, limits = c(llim, ulim), breaks = c(-2,0,2),guide = guide_colorbar(barwidth = 5, barheight = 1, 
                                                                                            ticks = TRUE, ticks.colour="black",ticks.linewidth=1/.pt,direction="horizontal",#"vertical",                                                                                            
                                                                                            label.theme = element_text(size = 8),
                                                                                            draw.ulim = TRUE, 
                                                                                            draw.llim = TRUE, 
                                                                                            draw.separator = TRUE,
                                                                                            
                                                                                            title.position = "top"
                                                                                            )        
            )

        pdf(file=paste0("write/Stat/Fig5_e_Pheatmap_ct4_DEP_limma_",Sys.Date(),".pdf"))
        print(p1)
        print(p2)
        dev.off()                                                                                            

    }
   
    # Desc: Pheatmap of sig.protein panel of 14 ct
    # Input: Obj_list$limma_ct14_df; PM proteins; Obj_list$impute_mean_df
    # Output: Pheatmap[Fig.6c]; Related table
    {
        rm(list=ls())
        Obj_list = readRDS(file="dataset_04.rds") 

        meta_df = Obj_list$meta_df %>% dplyr::select(CellType, CellType_Color) %>% unique()
        meta_df$CellType = factor(meta_df$CellType, levels = Obj_list$ct14_order)
        
        # filter sig. and PM proteins
        ct14_dep_df = Obj_list$limma_ct14_df %>% dplyr::filter(sum_pvals<0.05 & sum_log2fc>1)#5570
        row.names(ct14_dep_df) = ct14_dep_df$pids
        pm_df = Obj_list$anno_df %>% dplyr::filter(loci == "pm")#786                
        ct14_dep_pm_df = na.omit(ct14_dep_df[pm_df$pids,]) #475
                
        ct14_dep_pm_df$CellType = factor(ct14_dep_pm_df$CellType, levels = Obj_list$ct14_order)
        ct14_dep_pm_df = ct14_dep_pm_df %>% dplyr::arrange(CellType)

        write.csv(ct14_dep_pm_df, file=paste0("write/Stat/Fig6_c_Pheatmap_ct14_DEP_limma_mean_",Sys.Date(),".csv"))
        
        plot_df = Obj_list$impute_mean_df[ct14_dep_pm_df$pids, Obj_list$ct14_order] #475 * 14
                                  
        
        ## annotation
        annotation_col = data.frame(CellType = Obj_list$ct14_order)
        row.names(annotation_col) = Obj_list$ct14_order
        anno_colors = list(CellType = unique(Obj_list$ct14_color))
        names(anno_colors$CellType) = Obj_list$ct14_order
                
        # 使用apply函数对data.frame的每一行应用缩放函数  
        df_scaled_rows <- t(apply(plot_df, 1, scale))  
        
        # 计算99%分位数
        ulim <- quantile(df_scaled_rows, 0.99)
        llim <- quantile(df_scaled_rows, 0.01)

        # # 将超过99%分位数的值替换为99%分位数
        df_scaled_rows <- ifelse(df_scaled_rows > ulim, ulim, df_scaled_rows)
        df_scaled_rows <- ifelse(df_scaled_rows < llim, llim, df_scaled_rows)

        # 将结果转换回data.frame格式  
        df_scaled_rows <- as.data.frame(df_scaled_rows)
        names(df_scaled_rows) = names(plot_df)

        length(which(df_scaled_rows>0)) #62
        length(which(df_scaled_rows<0)) #71
        # 创建颜色映射
        my_palette <- c(
            c(colorRampPalette(c("navy", "white"))(29)),
            c(colorRampPalette(c("white", "red"))(37))
        )
        # 将0设置为断点
        my_breaks <- unique(c(seq(min(df_scaled_rows), 0, length.out = 29), seq(0, max(df_scaled_rows), length.out = 37)))

        p1=pheatmap::pheatmap(df_scaled_rows,cluster_rows = F,show_rownames = F,cluster_cols = F,
                scale = "none",border_color = NA,legend = T,cellwidth = 20,cellheight = 0.7,
                # color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                color = my_palette,breaks=my_breaks,#colorRampPalette(c("navy", "white", "red"))(100),
                annotation_col = annotation_col,show_colnames = F,legend_labels = F,annotation_legend = F,labels_col = NA,labels_row = NA,
                annotation_names_col = FALSE,annotation_colors = anno_colors,silent = TRUE
        )

        # Fig lengends           
        df <- data.frame(matrix(nrow = 10, ncol = 2))
        df[] <- rnorm(20)
        p2=ggplot(data = df, aes(x = X1, y = X2, fill = X2)) + geom_tile() + theme_void()+
            scale_fill_gradientn(colors=my_palette, limits = c(llim, ulim), breaks = c(-1.5,0,1.5,2.5),guide = guide_colorbar(barwidth = 5, barheight = 1, 
                                                                                            ticks = TRUE, ticks.colour="black",ticks.linewidth=1/.pt,direction="horizontal",#"vertical",                                                                                            
                                                                                            label.theme = element_text(size = 8),
                                                                                            draw.ulim = TRUE, 
                                                                                            draw.llim = TRUE, 
                                                                                            draw.separator = TRUE,
                                                                                            
                                                                                            title.position = "top"
                                                                                            )        
            )

        pdf(file=paste0("write/Stat/Fig6_c_Pheatmap_ct14_DEP_limma_mean_",Sys.Date(),".pdf"),height=8)
        print(p1)
        print(p2)
        dev.off()                                                                                            

    }

    # Desc: Line plot of sig. and surface protein panel of 14 ct
    # Input: Obj_list$limma_ct14_df; PM proteins; Obj_list$impute_mean_df
    # Output: Line plot[Fig.6d]; Related table
    # NOT FINISH
    {

    }



}

# -------------------------------------------------------------------------
#                          IV-Functional Enrichment analysis              #
# -------------------------------------------------------------------------
{
    # Desc: GO enrichment of cell-type enriched proteins for four major cell-lineages
    # Input: 
    # Output: 
    {
        rm(list=ls())        
        Obj_list = readRDS(file="dataset_04.rds") 

        meta_df = Obj_list$meta_df                
        data_df = read.delim(file="write/Stat/for_Fig5_f_GO-BP.csv",sep=",",header=T,row.names=1)

        data_df$pid = gsub(data_df$pids,pattern=";.*",replacement="")        
        tran_df = bitr(data_df$pids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        plot_df = data_df %>% merge(tran_df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)    
        plot_df$CellLineage = factor(plot_df$CellLineage, levels=Obj_list$ct4_order)    

        formula_res <- compareCluster(data = plot_df,ENTREZID~CellLineage,fun="enrichGO", 
                                    OrgDb = org.Mm.eg.db,
                                    ont = "BP", pAdjustMethod = "BH",
                                    pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                    
        )
        formula_res_cutoff = formula_res
        formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<=0.05,]

        write.csv(x = formula_res@compareClusterResult, file=paste0("write/Stat/Fig5_f_GOBP_",Sys.Date(),".csv"))  
        
        selected_GO = formula_res_cutoff@compareClusterResult
        selected_items = c("mitochondrial gene expression","translation","mitochondrial translation",                                                      
                           "extracellular matrix organization","extracellular structure organization","external encapsulating structure organization",
                           "B cell activation","leukocyte activation","B cell proliferation","inflammatory response",
                           "defense response to other organism","extracellular structure organization","inflammatory response","leukocyte activation"
                           
        )
        selected_GO = selected_GO %>% filter(Description %in% selected_items)
        formula_res_cutoff@compareClusterResult = selected_GO
        formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)

        p1=dotplot(formula_res_cutoff,color="pvalue",showCategory = 3,font.size=12,size="Count") 
        pdf(file=paste0("write/Stat/Fig5_f_GOBP_Selected_compareCluster_",Sys.Date(),".pdf"),width = 8,height = 6)
        print(p1)
        dev.off() 

        formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)
        p2=enrichplot::dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=14,color="log10P",size="count") + theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) #+
        # scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, 
        #                                                                         name = "YlOrRd")[3:6])(30))
        pdf(file=paste0("write/Stat/Fig5_f_GOBP_Selected_compareCluster2_",Sys.Date(),".pdf"),width = 8,height = 6)
        print(p2)
        dev.off()         
    }

}
