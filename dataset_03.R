#'@time 20240222
#'@author Yuan
#'@desc Bioinformastic analysis for spaital-proteomics after iPAC technology. 

library(dplyr)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
setwd("")

# -------------------------------------------------------------------------
# Part I, Quality Control and Pre-treat Quantified_Proteins files
{
    # Desc: Create .RDS and save files from Perseus output, including Quantified, Filter and Impute
    # Output: .RDS
    {
        rm(list=ls())
        Obj_list = list() # save RDS
        RDS_file_name = "obj_list_ct6_20240222.rds"

        # prepare annotation files
        anno_df = read.delim("spatial_afterNorm_Raw_6379.txt",header=T,check.names=FALSE) 
        anno_df = anno_df[,1:2] %>% unique()
        names(anno_df) = c("pids","genenames")
        anno_df$pid = gsub(anno_df$pids,pattern=";.*",replacement="")
        anno_df$genename = gsub(anno_df$genenames,pattern=";.*",replacement="")
        row.names(anno_df) = anno_df$pids

        # prepare meta files
        meta_df = read.delim("meta.txt",header=T,check.names=FALSE) 
        head(meta_df)

        # ct6 order
        ct6_order = c("Acinar", "PanIN", "PDAC", "CAF", "IT", "LN")
        ct3_order = c("Acinar", "PanIN", "PDAC")
        ct2_order = c("IT", "LN")

        # prepare Quantified file
        quantified_df = read.delim("Quantified_ct6_6377.txt",header=T,check.names=FALSE) #6377 * 18    
        quantified_df = quantified_df[,meta_df$ShortSample]#re-order
        names(quantified_df) = meta_df$SampleID #re-name

        # prepare Filter file
        filter_df = read.delim("Filter_ct6_5844.txt",header=T,check.names=FALSE) #5844 * 18    
        row.names(filter_df) = filter_df$Protein
        filter_df = filter_df[,meta_df$ShortSample]
        names(filter_df) = meta_df$SampleID

        # prepare Impute file
        impute_df = read.delim("Impute_ct6_5844.txt",header=T,check.names=FALSE) #5844 * 18    
        row.names(impute_df) = impute_df$Protein
        impute_df = impute_df[,meta_df$ShortSample]
        names(impute_df) = meta_df$SampleID

        # prepare Impute_Median file
        impute_median_df = read.delim("Impute_median_ct6_5844.txt",header=T,check.names=FALSE) #5844 * 18    
        row.names(impute_median_df) = impute_median_df$Protein
        impute_median_df = impute_median_df[,ct6_order]

        Obj_list$quantified_df = quantified_df
        Obj_list$filter_df = filter_df
        Obj_list$impute_df = impute_df
        Obj_list$impute_median_df = impute_median_df
        Obj_list$meta_df = meta_df
        Obj_list$anno_df = anno_df
        Obj_list$ct6_order = ct6_order
        Obj_list$ct3_order = ct3_order
        Obj_list$ct2_order = ct2_order

        saveRDS(Obj_list, file=RDS_file_name)

        # for tangram
        Obj_list = readRDS(file = "obj_list_ct6_20240222.rds")
        write.csv(Obj_list$meta_df,file="/aaa/zihanwu/yyyli2/projectx_xu_NC2024/raw/Tangram_deconvolution/spatial6_filter_df_5844_meta.csv")
        filter_df = Obj_list$filter_df[Obj_list$meta_df$SampleID]
        row.names(filter_df) = gsub(row.names(filter_df), pattern=";.*", replacement="")
        write.csv(filter_df,file="/aaa/zihanwu/yyyli2/projectx_xu_NC2024/raw/Tangram_deconvolution/spatial6_filter_df_5844.csv")
    }

    # Desc: for PCA analysis in Spectronaut
    # Input: Obj_list$quantified_df
    # for Fig.4b,f
    # Anno: may be merged 
    {
        rm(list=ls()) 
        Obj_list = readRDS(file="obj_list_ct6_20240222.rds")

        quantified_df = Obj_list$quantified_df
        quantified_df$pids = row.names(quantified_df)
        output_df = merge(quantified_df, Obj_list$anno_df,by="pids",all=F)
        dim(output_df) #6377 * 22

        write.csv(output_df,
                    file = paste0("write/Quantified/Quantified_for_PCA_",Sys.Date(),".csv"),row.names = T)    

    }

    # Desc: ANOVA analysis for cell-type enrichment and functional analysis, including Acinar PanIN and PDAC
    # Input: Obj_list$impute_df; Obj_list$impute_median_df
    # Output: Obj_list$aov_ct3_df
    # for Fig.4c,d,e
    {
        rm(list=ls())
        Obj_list = readRDS(file="obj_list_ct6_20240222.rds") 

        meta_df = Obj_list$meta_df %>% filter(Group %in% Obj_list$ct3_order)# 12 * 4
        meta_df$celltype = meta_df$Group
        meta_df$celltype = factor(meta_df$celltype, levels = Obj_list$ct3_order)

        protein_df = Obj_list$impute_df[,meta_df$SampleID] #5844 * 9
        protein_median_df = Obj_list$impute_median_df[,Obj_list$ct3_order]
        
        protein_num = dim(protein_df)[1]
        sample_num = dim(protein_df)[2]

        ## DEP with median  
        dep_df = data.frame()
        for(i in 1:3){
            temp_df = protein_median_df[,i] - protein_median_df[,1:3]
            temp_pid = names(which(rowSums(temp_df>log2(1)) == 2)) # must large than the rest groups
            temp_df2 = data.frame(pids=temp_pid)
            temp_df2$celltype = names(protein_median_df)[i]
            temp_df2$log2FC = apply(temp_df[temp_pid,-i],1,median)        

            dep_df = rbind(dep_df,temp_df2)            
        }    


        Pvalue = c()
        for(i in 1:protein_num){
            pid_df = protein_df[i,] %>% t() %>% as.data.frame()
            names(pid_df) = "pid"
            pid_df$SampleId = row.names(pid_df)
            pid_df$celltype = meta_df$celltype
            temp_pvalue = summary(aov(data=pid_df,pid~celltype))[[1]]$`Pr(>F)`[1]
            Pvalue = c(Pvalue,temp_pvalue)
            
        }

        protein_df$pvalue = Pvalue
        protein_df$fdr = p.adjust(Pvalue,method = "BH")
        protein_df$pids = row.names(protein_df)

        length(which(protein_df$fdr <0.05)) #698
        length(which(protein_df$pvalue <0.05)) #1619
        
        
        # output
        output_df = merge(protein_df,dep_df,by="pids",all=F)
        Obj_list$aov_ct3_df = output_df
        saveRDS(Obj_list, file="obj_list_ct6_20240222.rds")


        output_df = merge(Obj_list$anno_df,output_df,by="pids",all=F)
        write.csv(output_df,
                    file = paste0("write/Stat/ANOVA_for_enrich_",Sys.Date(),".csv"),row.names = T)   


    }    

    # Desc: Significant analysis of cell-type sig. diff. proteins, including IT and LN
    # Input: Obj_list$quantified_df
    # Output: Obj_list$ttest_ct2_df
    # for Fig.4g,h
    {
        rm(list=ls())
        Obj_list = readRDS(file="obj_list_ct6_20240222.rds") 

        meta_df = Obj_list$meta_df %>% filter(Group %in% Obj_list$ct2_order)# 6 * 4
        meta_df$celltype = meta_df$Group
        meta_df$celltype = factor(meta_df$celltype, levels = Obj_list$ct2_order)

        protein_df = Obj_list$impute_df[,meta_df$SampleID] #5844 * 9
        # protein_median_df = Obj_list$impute_median_df[,Obj_list$ct2_order]

        ct_A_ind = which(meta_df$celltype == Obj_list$ct2_order[1])
        ct_B_ind = which(meta_df$celltype == Obj_list$ct2_order[2])     
        Pvalue = apply(protein_df,MARGIN = 1,FUN = function(x){
            return(t.test(x[ct_A_ind],x[ct_B_ind])$p.value)
            })

        log2FC = apply(protein_df,MARGIN = 1,FUN = function(x){
            y = median(x[ct_A_ind]) - median(x[ct_B_ind])
            return(y)
        })        

        protein_df$pvalue = Pvalue
        protein_df$fdr = p.adjust(Pvalue,method = "BH")
        protein_df$log2FC = log2FC

        protein_df$pids = row.names(protein_df)

        length(which(protein_df$fdr <0.05)) #0
        length(which(protein_df$pvalue <0.05)) #799
        
        # output
        output_df = protein_df
        Obj_list$ttest_ct2_df = output_df
        saveRDS(Obj_list, file="obj_list_ct6_20240222.rds")


        output_df = merge(Obj_list$anno_df,output_df,by="pids",all=F)
        write.csv(output_df,
                    file = paste0("write/Stat/ttest_for_sig_",Sys.Date(),".csv"),row.names = F)   
    
    
    }
}

# -------------------------------------------------------------------------
# Part II, Bioinformastic analysis of all 6 cell types
{
    # Desc: Box plot for quantified protein number
    # Input: Obj_list$quantified_df
    # Fig.4a
    {
        rm(list=ls()) 
        Obj_list = readRDS(file="obj_list_ct6_20240222.rds")
        valid_num = apply(Obj_list$quantified,2,FUN = function(x){
            length(which(as.numeric(x)!=0))
        })

        plot_df = data.frame(protein_count = valid_num, group=Obj_list$meta$Group)
        plot_df$group = factor(plot_df$group, levels = Obj_list$ct6_order)

        p1=ggplot(data = plot_df,mapping=aes(x=group,y=protein_count,fill=group),color="black")+
            geom_boxplot(lwd=0.5,outlier.size = 0.5) +
            stat_boxplot(geom = "errorbar",#color="black",
                        lwd=0.5, 
                        width=0.2)+
            geom_jitter(size=1, alpha=10,width = 0.2) +
            theme_bw()+
            theme(panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                axis.title = element_text(face="bold"),
                axis.text.x = element_text(    angle = 45,
                                                hjust = 0.5,
                                                vjust = 0.5),
                # axis.ticks.x  = element_blank(),
                axis.text.y = element_text(colour="black", size = 9),
                axis.line = element_line(size=0.5, colour = "black"))+
            scale_y_continuous(limits=c(0,6000),breaks=seq(0,6000,1000),label=seq(0,6000,1000),expand = c(0,1))+
            labs(x="",y="Quantified Protein Number")
        pdf(paste0("write/Quantified/Quantified_Box_",Sys.Date(),".pdf"),width = 4,height = 3)
        print(p1)
        dev.off()

        ## save data
        plot_df$SampleID = row.names(plot_df)
        output_df = merge(Obj_list$meta, plot_df, by="SampleID",all=F)
        output_df$group = factor(output_df$group, levels=Obj_list$ct6_order)
        output_df = output_df %>% arrange(group)
        head(output_df)
        write.csv(output_df,
                    file = paste0("write/Quantified/Quantified_Box_",Sys.Date(),".csv"),row.names = T)
    }
            
}

# -------------------------------------------------------------------------
# Part III, Bioinformastic analysis of PDAC development related 3 cell types, including Acinar PanIN and PDAC
{
    # Desc: Pheatmap of selected proteins, including Acinar PanIN and PDAC
    # Input: Obj_list$aov_ct3_df
    # Output: Pheatmap[Fig.4c]; Related table
    {
        rm(list=ls())
        Obj_list = readRDS(file="obj_list_ct6_20240222.rds") 

        meta_df = Obj_list$meta_df %>% filter(Group %in% Obj_list$ct3_order)# 12 * 4
        meta_df$celltype = meta_df$Group
        meta_df$celltype = factor(meta_df$celltype, levels = Obj_list$ct3_order)

        plot_df =  protein_df[,meta_df$ShortSample] # use impute data
        row.names(plot_df) = protein_df$Protein   
        plot_df = plot_df[dep_df$pid,] # select sig pid
        
        bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
        p1=pheatmap::pheatmap(plot_df,scale="row",cluster_rows = F,cluster_cols = T,show_colnames=T,
                    show_rownames = F,#border_color = "white",
                    #   legend_breaks=seq(-1,1,1),#kmeans_k = 20,
                    cellwidth = 10,cellheight = 0.3,silent = T,
                    #   annotation_col = plot_anno_df,#annotation_colors = color.list,
                    #   color = c(colorRampPalette(colors = c("blue","navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3","red"))(length(bk)/2))
                    # color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2))
                    color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
        )
        
        pdf(paste0("iPAC_KPC_application/Pheatmap_",Sys.Date(),".pdf"),width=8,height = 10)
        print(p1)
        dev.off()

        # output_df = merge(dep_df, protein_df,by.y="Protein",by.x="pid",all=F)
        # write.table(output_df, file = paste0("iPAC_KPC_application/ANOVA_for_pheatmap_",Sys.Date(),".txt"),sep="\t",row.names=FALSE)

    }

    # Desc: GO enrichment of cell-type enriched proteins, including Acinar PanIN and PDAC
    # Input: Obj_list$aov_ct3_df
    # Output: Dot plot[Fig.4d]; Related table
    {
        rm(list=ls())


        Obj.list = readRDS(file = "write/KPC_Obj_20230505.rdata")
        dep.df = Obj.list$PDAC_stage %>% filter(Pvalue<0.05 & abs(log2FC)>=log2(2))
        dep.df$ENTREZID = Obj.list$anno[dep.df$pid,"GeneID"]
        # tran.df = bitr(dep.df$pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        # dep.df = dep.df %>% merge(tran.df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)

        # bg_pid = row.names(Obj.list$impute)
        # bg_pid.df = bitr(bg_pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

        formula_res <- compareCluster(data = dep.df,ENTREZID~region,fun="enrichGO", 
                                    OrgDb = org.Mm.eg.db,
                                    ont = "BP", pAdjustMethod = "BH",
                                    pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE,
                                    universe = Obj.list$anno$GeneID
        )
        formula_res_cutoff = formula_res
        formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<=0.05,]

        write.csv(x = formula_res_cutoff@compareClusterResult, file=paste0("write/02_PDAC/pvalue_sig_GOBP_",Sys.Date(),".csv"))  

        # 
        # selected_GO = formula_res_cutoff@compareClusterResult
        # selected_items = c("ribonucleoprotein complex biogenesis",
        #                    "ncRNA metabolic process",
        #                    "mRNA processing",
        #                    "mononuclear cell differentiation",
        #                    "lymphocyte differentiation",
        #                    "myeloid leukocyte activation",
        #                    "extracellular matrix organization",
        #                    "extracellular structure organization",
        #                    "external encapsulating structure organization"
        # )
        # selected_GO = selected_GO %>% filter(Description %in% selected_items)
        formula_res_cutoff@compareClusterResult = selected_GO
        formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)
        p2=dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=14,color="log10P",size="count") + theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
        scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, 
                                                                                name = "YlOrRd")[3:6])(30))
        pdf(file=paste0("write/02_PDAC/GOBP_Selected_compareCluster_",Sys.Date(),".pdf"),width = 8,height = 6)
        print(p2)
        dev.off() 
    }

    # Desc: Sig analysis of cell-type enriched proteins, including Acinar PanIN and PDAC
    # Input: Obj_list$filter
    # Output: Box plot[Ext Data Fig.9]; Related table
    {
        rm(list=ls())
        
        sig_input.df = input.df %>% filter(input.df$Pvalue<0.05)
        box_spread.df = sig_input.df[,-c(10,11)] %>% tidyr::gather(key="samples",value="intensity",-entry)
        box_spread.df$group = gsub(box_spread.df$samples,pattern = "_.*$",replacement = "")
        box_spread.df$group = factor(box_spread.df$group, levels = unique(box_spread.df$group))
        
        
        list_merge.list= list()
        
        comparison_list = list(c("Acinar","PanIN"),c("Acinar","PDAC"),c("PanIN","PDAC"))
        
        for (protein in unique(sig_input.df$entry)){
        # print(protein)
        
        temp_box.df = box_spread.df %>% dplyr::filter(entry == protein)
        
        
        p1=ggplot(data = temp_box.df,mapping=aes(x=group,y=intensity,fill=group),color="black")+
            geom_boxplot(lwd=0.5,outlier.size = 0.5) +
            stat_boxplot(geom = "errorbar",#color="black",
                        lwd=0.5, 
                        width=0.2)+
            geom_jitter(size=1, alpha=10,width = 0.2) +
            theme_classic()+
            theme(panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                # text=element_text(family = "Tahoma", face ="bold"),
                axis.title = element_text(face="bold"),
                # axis.text.x = element_text(    angle = 90,
                #                                 hjust = 1,
                #                                 vjust = 0.5),
                # axis.ticks.x  = element_blank(),
                axis.text.y = element_text(colour="black", size = 9),
                axis.text.x = element_text(colour="black", size = 9),
                axis.line = element_line(size=0.5, colour = "black"))+
            # scale_y_continuous(limits=c(0,20),breaks=seq(0,20,2),label=seq(0,20,2),expand = c(0,1))+
            labs(x="",y="log2 LFQ intensity",title = protein)+
            stat_compare_means(comparisons = comparison_list, label="p.signif", paired = T,method="t.test")
        
        list_merge.list[[protein]] = p1
        }
        
        merge.plot = ggarrange(plotlist = list_merge.list,ncol = 1,nrow = 1)
        pdf(paste0("write/Quantified/Merged_PDAC_Sig_Quantified_Box_",Sys.Date(),".pdf"),width = 4,height = 3,onefile = T)
        print(merge.plot)
        dev.off()
        
        
        ## save data
        write.csv(box_spread.df,
                file = paste0("write/Quantified/Merged_PDAC_Sig_Quantified_Box_",Sys.Date(),".csv"),row.names = T)
    
    
    }

}

# -------------------------------------------------------------------------
# Part IV, Bioinformastic analysis of PDAC development related 3 cell types, including IT and LN
{
    # Desc: Pheatmap of selected proteins, including IT and LN
    # Input: Obj_list$ttest_ct2_df
    # Output: Pheatmap[Fig.4g]; Related table
    {
        rm(list=ls())
        Obj_list = readRDS(file="obj_list_ct6_20240222.rds") 

        meta_df = Obj_list$meta_df %>% filter(Group %in% Obj_list$ct3_order)# 12 * 4
        meta_df$celltype = meta_df$Group
        meta_df$celltype = factor(meta_df$celltype, levels = Obj_list$ct3_order)

        plot_df =  protein_df[,meta_df$ShortSample] # use impute data
        row.names(plot_df) = protein_df$Protein   
        plot_df = plot_df[dep_df$pid,] # select sig pid
        
        bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
        p1=pheatmap::pheatmap(plot_df,scale="row",cluster_rows = F,cluster_cols = T,show_colnames=T,
                    show_rownames = F,#border_color = "white",
                    #   legend_breaks=seq(-1,1,1),#kmeans_k = 20,
                    cellwidth = 10,cellheight = 0.3,silent = T,
                    #   annotation_col = plot_anno_df,#annotation_colors = color.list,
                    #   color = c(colorRampPalette(colors = c("blue","navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3","red"))(length(bk)/2))
                    # color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2))
                    color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))
        )
        
        pdf(paste0("iPAC_KPC_application/Pheatmap_",Sys.Date(),".pdf"),width=8,height = 10)
        print(p1)
        dev.off()

        # output_df = merge(dep_df, protein_df,by.y="Protein",by.x="pid",all=F)
        # write.table(output_df, file = paste0("iPAC_KPC_application/ANOVA_for_pheatmap_",Sys.Date(),".txt"),sep="\t",row.names=FALSE)

    }

    # Desc: GO enrichment of cell-type enriched proteins, including IT and LN
    # Input: Obj_list$ttest_ct2_df
    # Output: Dot plot[Fig.4h]; Related table
    {
        rm(list=ls())
        library(enrichplot)
        library(clusterProfiler)
        library(org.Mm.eg.db)

        Obj.list = readRDS(file = "write/KPC_Obj_20230505.rdata")
        dep.df = Obj.list$PDAC_stage %>% filter(Pvalue<0.05 & abs(log2FC)>=log2(2))
        dep.df$ENTREZID = Obj.list$anno[dep.df$pid,"GeneID"]
        # tran.df = bitr(dep.df$pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        # dep.df = dep.df %>% merge(tran.df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)

        # bg_pid = row.names(Obj.list$impute)
        # bg_pid.df = bitr(bg_pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

        formula_res <- compareCluster(data = dep.df,ENTREZID~region,fun="enrichGO", 
                                    OrgDb = org.Mm.eg.db,
                                    ont = "BP", pAdjustMethod = "BH",
                                    pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE,
                                    universe = Obj.list$anno$GeneID
        )
        formula_res_cutoff = formula_res
        formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<=0.05,]

        write.csv(x = formula_res_cutoff@compareClusterResult, file=paste0("write/02_PDAC/pvalue_sig_GOBP_",Sys.Date(),".csv"))  

        # 
        # selected_GO = formula_res_cutoff@compareClusterResult
        # selected_items = c("ribonucleoprotein complex biogenesis",
        #                    "ncRNA metabolic process",
        #                    "mRNA processing",
        #                    "mononuclear cell differentiation",
        #                    "lymphocyte differentiation",
        #                    "myeloid leukocyte activation",
        #                    "extracellular matrix organization",
        #                    "extracellular structure organization",
        #                    "external encapsulating structure organization"
        # )
        # selected_GO = selected_GO %>% filter(Description %in% selected_items)
        formula_res_cutoff@compareClusterResult = selected_GO
        formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)
        p2=dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=14,color="log10P",size="count") + theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
        scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, 
                                                                                name = "YlOrRd")[3:6])(30))
        pdf(file=paste0("write/02_PDAC/GOBP_Selected_compareCluster_",Sys.Date(),".pdf"),width = 8,height = 6)
        print(p2)
        dev.off() 
    }

}
