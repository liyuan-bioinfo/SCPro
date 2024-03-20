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
setwd("/aaa/zihanwu/yyyli2/projectx_xu_NC2024/raw/iPAC_KPC_ct8_valid")

# -------------------------------------------------------------------------
#                            I-Basic pre-treat                            #
# -------------------------------------------------------------------------
{
    # Desc: Basic check ct8 files
    # Input: ProteinGroups_Quantified file from MS-search software
    # Output: Quantified_Proteins.txt    
    {        
        rm(list=ls())
        Obj_list = readRDS("Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")

        names(Obj_list)

        data_df = read.delim("KLRG1_allsubtypes_combined_protein(1).tsv")

        head(Obj_list$identified)
        write.csv(Obj_list$meta, file="meta.csv",row.names=FALSE)
    }

    # Desc: for RDS files
    {
        rm(list=ls())
        protein_data_name = "KLRG1_allsubtypes_combined_protein(1).tsv"
        meta_data_name = "meta.csv"    
        
        protein_df = read.delim(protein_data_name,header=T,check.names=FALSE,sep="\t") #5325 * 232
        meta_df = read.delim(meta_data_name,header=T,sep=",") #24 * 4
       
        
        ## quantified proteins    
        quantified_protein_df = protein_df[,meta_df$Sample] #5153 * 24
        row.names(quantified_protein_df) = protein_df$`Protein ID`
        names(quantified_protein_df) = meta_df$SampleID
        
        quantified_protein_df <- replace(quantified_protein_df, is.na(quantified_protein_df), 0)

        temp_sum = rowSums(quantified_protein_df, na.rm = TRUE)
        quantified_protein_df = quantified_protein_df[which(temp_sum != 0 ), ] #4077 * 24

        ## save Quantified table for Perseus    
        write.table(quantified_protein_df, file = "Quantified_ct8_4077.txt",sep="\t")    

        ## identified proteins
        meta_df$Sample = gsub(meta_df$Sample, pattern = "MaxLFQ ",replacement="")
        identified_protein_df = protein_df[,meta_df$Sample] #5153 * 24
        row.names(identified_protein_df) = protein_df$`Protein ID`
        names(identified_protein_df) = meta_df$SampleID
        
        identified_protein_df <- replace(identified_protein_df, is.na(identified_protein_df), 0)

        temp_sum = rowSums(identified_protein_df, na.rm = TRUE)
        identified_protein_df = identified_protein_df[which(temp_sum != 0 ), ] #5029 * 24
        write.table(identified_protein_df, file = "Identified_ct8_5029.txt",sep="\t")    
        
        
    }    
    # Desc: Create .RDS and save files from Perseus output, including Quantified, Filter and Impute
    {
        rm(list=ls())
        Obj_list = list() # save RDS
        RDS_file_name = "obj_list_ct8_20240316.rds"

        # prepare annotation files
        anno_df = read.delim("KLRG1_allsubtypes_combined_protein(1).tsv",header=T,check.names=FALSE,sep="\t") 
        anno_df = anno_df[,c("Protein ID","Gene")] %>% unique()
        names(anno_df) = c("pids","genenames")
        anno_df$pid = gsub(anno_df$pids,pattern=";.*",replacement="")
        anno_df$genename = gsub(anno_df$genenames,pattern=";.*",replacement="")
        row.names(anno_df) = anno_df$pids #5153

        # prepare meta files
        # ct8 order
        ct8_order = c("KPC_425KN","KPC_425KP","KPC_425N","KPC_425P","KPC_4KN","KPC_4KP","KPC_8KN","KPC_8KP")
        group_order = c("KPC_425K","KPC_425","KPC_4K","KPC_8K")
        ct8_color = c("#08306B","#08519C","#2171B5","#4292C6","#FD8D3C","#F16913","#ADDD8E","#78C679")
        meta_df = read.delim("meta.csv",header=T,check.names=FALSE,sep=",") 
        meta_df$CellType_color = factor(meta_df$CellType_color, levels=ct8_color)
        meta_df$CellType = factor(meta_df$CellType, levels=ct8_order)
        meta_df$Group = factor(meta_df$Group, levels=group_order)
        meta_df = meta_df %>% arrange()
        head(meta_df)
                
        # prepare Identified file
        identified_df = read.delim("Identified_ct8_5029.txt",header=T,check.names=FALSE) #5029 * 24
        identified_df = identified_df[,meta_df$SampleID]#re-order
        
        # prepare Quantified file
        quantified_df = read.delim("Quantified_ct8_4077.txt",header=T,check.names=FALSE) #4077 * 24
        quantified_df = quantified_df[,meta_df$SampleID]#re-order        

        # prepare Filter file
        filter_df = read.delim("Filter_ct8_3824.txt",header=T,check.names=FALSE) #3824 * 24
        row.names(filter_df) = filter_df$pid
        filter_df = filter_df[,meta_df$SampleID]        

        # prepare Impute file
        impute_df = read.delim("Impute_ct8_3824.txt",header=T,check.names=FALSE) #3824 * 24
        row.names(impute_df) = impute_df$pid
        impute_df = impute_df[,meta_df$SampleID]
        
        # set color of ct8
        Obj_list$identified_df = identified_df
        Obj_list$quantified_df = quantified_df
        Obj_list$filter_df = filter_df
        Obj_list$impute_df = impute_df
        
        Obj_list$meta_df = meta_df
        Obj_list$anno_df = anno_df
        Obj_list$ct8_order = ct8_order        
        Obj_list$ct8_color = ct8_color


        saveRDS(Obj_list, file=RDS_file_name)        

        # for tangram
        Obj_list = readRDS(file = "obj_list_ct8_20240316.rds")
        write.csv(Obj_list$meta_df,file="/aaa/zihanwu/yyyli2/projectx_xu_NC2024/raw/Tangram_deconvolution/ct8_filter_df_3824_meta.csv")
        write.csv(Obj_list$filter_df[Obj_list$meta_df$SampleID],file="/aaa/zihanwu/yyyli2/projectx_xu_NC2024/raw/Tangram_deconvolution/ct8_filter_df_3824.csv")
    }

    # Desc: Stat analysis
    # Student's test for two-tail sample statistical analysis
    {
        rm(list=ls())
        obj_list = readRDS(file = "obj_list_ct8_20240316.rds")
        meta_df = obj_list$meta_df %>% filter(Group == "KPC_425K")
        meta_df$CellType = factor(meta_df$CellType, levels=c("KPC_425KP", "KPC_425KN"))
        meta_df = meta_df %>% arrange(CellType)

        data_df = obj_list$impute_df[,meta_df$SampleID]
        protein_num = dim(data_df)[1]                
        
        Pvalue = c()
        log2FC = c()
        for(i in 1:protein_num){
            temp_df = data_df[i,] %>% t() %>% as.data.frame()            
            names(temp_df) = "intensity"
            temp_df$group = meta_df$CellType

            # temp_df <- temp_df[!(temp_df$intensity == 0), ]
            model = t.test(intensity ~ group, data=temp_df)
            temp_pvalue = model$p.value
            temp_log2fc = model$estimate[[1]] - model$estimate[[2]]

            Pvalue = c(Pvalue,temp_pvalue)
            log2FC = c(log2FC,temp_log2fc) # Group A - Group B
            
        }
        data_df$pvalue = Pvalue
        data_df$fdr = p.adjust(Pvalue,method = "BH")
        data_df$log2fc = log2FC
        
        data_df$pid = row.names(data_df)
        data_df$genename = obj_list$anno_df[data_df$pid,"genename"]
        data_df$sig = ""

        data_df[which(data_df$pvalue <0.05 & data_df$log2fc>1),"sig"] = "up"
        data_df[which(data_df$pvalue <0.05 & data_df$log2fc<=-1),"sig"] = "down"
        table(data_df$sig)
        
        obj_list$dep_df = data_df
        obj_list$ct2_order = c("KPC_425KP", "KPC_425KN")
        saveRDS(obj_list,file="obj_list_ct8_20240316.rds")
        write.csv(x=data_df, file="write/Stat/ct8_t_test.csv",row.names = F)

    }

}

# -------------------------------------------------------------------------
#                            II-Quality Control                           #
# -------------------------------------------------------------------------
{
    # Desc: Bar plot for quantified protein number
    # Input: Obj_list$quantified_df and identified_df
    # Output: Bar plot [Extended Data Fig. 12]
    {
        rm(list=ls())
        Obj_list = readRDS(file = "obj_list_ct8_20240316.rds")
        temp.df = Obj_list$identified_df %>% t() %>% as.data.frame()
        
        meta.df = Obj_list$meta_df %>% dplyr::select(SampleID,CellType) %>% unique()
        # row.names(plot.data) = plot.data$group
        plot.data = meta.df %>% dplyr::select(CellType) %>% unique()
        row.names(plot.data) = plot.data$CellType
        
        temp_sum = c()
        for(i in Obj_list$ct8_order){
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
        plot.data2 = meta.df %>% dplyr::select(CellType) %>% unique()
        row.names(plot.data2) = plot.data2$CellType
        temp.df2 = Obj_list$quantified_df %>% t() %>% as.data.frame()
        
        temp_sum = c()
        for(i in Obj_list$ct8_order){
            temp_index = which(meta.df$CellType==i)  
            tt = temp.df2[temp_index,]
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
        plot.data.bar$CellType = factor(x = plot.data.bar$CellType,levels = c(Obj_list$ct8_order,"Total"))
        
        p1 <- ggplot(plot.data.bar, aes(x=CellType,y=protein_sum,fill=protein_loca)) + 
        geom_bar(stat = "identity",width=0.8)+
        scale_fill_grey(start = 0.75, end = 0.5) +
        theme_classic()+
        theme(        
            plot.title = element_text(size=16,hjust = 0.5),
            axis.title.y.left  = element_text(size=16),
            axis.text = element_text(size=12),
            axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1,color = Obj_list$ct8_color),
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
                    file = paste0("write/Quantified/sFig12_c_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".csv"),
                    fileEncoding = "UTF-8",row.names = F,quote = TRUE)
        
        pdf(paste0("write/Quantified/sFig12_c_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".pdf"),
                width=6,height = 4)
        print(p1)
        dev.off()
    }

    # Desc: Quantified analysis of 8 valid celltypes
    # Input: Obj_list$quantified
    # Output: PCA plot[Fig.6f]
    {  
        rm(list=ls())
        Obj_list = readRDS(file = "obj_list_ct8_20240316.rds")
        meta_df = Obj_list$meta_df
        
        data_df = Obj_list$quantified_df[,meta_df$SampleID]#select and sort
        data_df = log2(data_df+1) %>% t() %>% as.data.frame()
                
        pca_result <- prcomp(data_df,scale=T,center = T)        
        plot_df = as.data.frame(pca_result$x)

        summ1 <- summary(pca_result)
        xlab1 <- paste0("PC1 (",round(summ1$importance[2,1]*100,2),"%)")
        ylab1 <- paste0("PC2 (",round(summ1$importance[2,2]*100,2),"%)")
                
        plot_df = cbind(plot_df, meta_df)
        plot_df$CellType = factor(plot_df$CellType, levels = Obj_list$ct8_order)        
        
        p1=ggplot(data = plot_df,aes(x = PC1,y = PC2,color = CellType,
                                shape=Group
                                ))+

            geom_point(size = 3.5)+
            labs(x = xlab1,y = ylab1,color = "",title = "",shape="")+
            guides(fill = "none")+
            theme_bw()+            
            scale_colour_manual(values = Obj_list$ct8_color)+
            theme(plot.background = element_blank(),legend.background = element_blank(),
                panel.background = element_blank(),panel.grid = element_blank(),
                axis.text = element_text(size = 12),axis.title = element_text(size = 16),
                legend.text = element_text(size = 16),
                plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))

        # save plot        
        pdf(file=paste0("write/Quantified/Fig6_f_PCA_ct8_Quantified_",Sys.Date(),".pdf"),
            width = 7,height=5)
        print(p1)
        dev.off()
        
        write.csv(x = plot_df,file=paste0("write/Quantified/Fig6_f_PCA_ct8_Quantified_",Sys.Date(),".csv"))
    }

}

# -------------------------------------------------------------------------
#                            III-Stat Analysis                            #
# -------------------------------------------------------------------------
{
    # Desc: Volcano plot for CD4+CD25Klrg1- vs CD4+CD25+Klrg1+
    # Input: Obj_list$dep_df
    # Output: Volcano Plot[Fig.6g]; Related table
    {
        rm(list=ls())
        obj_list = readRDS(file="obj_list_ct8_20240316.rds") 
        meta_df = obj_list$meta_df        
        dep_df = obj_list$dep_df

        log10_P_Value_cutoff <- -log10(0.05)
        logFC_cutoff <- 1
        dep_df$`-log10_P.Value` = -log10(dep_df$pvalue)

        p1 <- ggplot(data = dep_df,aes(x = log2fc,y = `-log10_P.Value`))+
            geom_point(data = subset(dep_df,sig==""),size=0.5,
                    col = 'gray',alpha = 0.4)+
            geom_point(data = subset(dep_df,sig=="up"),size=0.5,
                    col = 'firebrick3',alpha = 1)+#down
            geom_point(data = subset(dep_df,sig=="down"),size=0.5,
                    col = 'blue',alpha = 1)+#up
            # geom_point(data = subset(temp.df.select,!is.na(gene_plot) ),
            #            aes(size = abs(logFC)),col="black",alpha = 0.4)+
            
            theme_bw()+
            labs(x='log2(Fold change)',y='-log10(p-value)')+
            theme(legend.title = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.background = element_blank(),
                plot.title = element_text(hjust=0.5,size = 14),
                text = element_text(size=12),
                legend.position = 'none',
                axis.line = element_line(colour = "black"))+
            labs(subtitle=paste0(length(which(dep_df$sig=="down")),"_",length(which(dep_df$sig=="up"))))+
        
        geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),lty = 3,col = 'black',lwd = 0.4)+
        geom_hline(yintercept = log10_P_Value_cutoff,lty = 3,col = 'black',lwd = 0.4)+
        
        ggrepel::geom_text_repel(data = subset(dep_df,genename%in%c("Cd4", "Foxp3","Ikzf2","Casp1","Kdelr2","Map4k1","Klrg1")),aes(label = genename),size = 3,
                max.overlaps = 20,min.segment.length = 0,seed = 123)
        #   geom_text_repel(data = temp.df.select,
        #                   max.overlaps = 20,min.segment.length = 0,seed = 123,
        #                   aes(label = label_volcano),size = 3,col = 'black')
        
                                                    
        pdf(file=paste0("write/Stat/ct8_t_test_VolcanoPlot_",Sys.Date(),".pdf"),width = 4,height = 4)
        print(p1)
        dev.off()                                                                                                

    }

}

# -------------------------------------------------------------------------
#                          IV-Functional Enrichment analysis              #
# -------------------------------------------------------------------------
{

   #GSEA分析
    {
        rm(list=ls())        
        obj_list = readRDS(file = "obj_list_ct8_20240316.rds")
        
        meta_df = obj_list$meta_df        

        dep_df = obj_list$dep_df
        dep_df$`-logp` = -log10(dep_df$pvalue)
        dep_df$rank = dep_df$log2fc * (dep_df$`-logp`)
        dep_df = dep_df %>% arrange(desc(rank))
        
        temp_geneList = dep_df$rank
        names(temp_geneList) = dep_df$pid
        re = clusterProfiler::gseGO(geneList =temp_geneList,OrgDb = org.Mm.eg.db,
                                    keyType = "UNIPROT",minGSSize = 3,maxGSSize = 2000,
                                    by = "fgsea",pvalueCutoff = 1,seed=123
        )
        
        saveRDS(re, paste0("write/Functional/ct8_DEP_GSEA_",Sys.Date(),".rds"))
        write.csv(x = re@result, file=paste0("write/Functional/Fig6h_ct8_DEP_GSEA_",Sys.Date(),".csv"))  


        # PLOT
        # re = readRDS("output_immune/GSEA.rds")
        selected_items = c(
            "myeloid leukocyte activation","myeloid leukocyte cytokine production"
            )
            
        selected.re = re
        selected.re@result = selected.re@result %>% filter(Description %in% selected_items)
        selected.re@result$Description = factor(selected.re@result$Description,levels = selected_items)
        selected.re@result = selected.re@result %>% arrange(Description)
                
        p1=gseaplot2(selected.re, "GO:0002274", base_size=14,pvalue_table = T,col = "firebrick")# + theme_bw()
        p2=gseaplot2(selected.re, "GO:0061082", base_size=14,pvalue_table = T,col = "firebrick")# + theme_bw()        
            
        pdf(file=paste0("write/Functional/Fig6h_ct8_DEP_GSEA_",Sys.Date(),".pdf"),width = 8,height = 6)
        print(p1)
        print(p2)
        dev.off() 
    
    }
}
