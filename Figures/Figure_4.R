library(readxl)
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

setwd("")
# Figure_4c, heat map of Sig. proteins for acinar, PanIN and PDAC regions
{
    rm(list=ls())
    Obj.list = readRDS(file="sp_dataset_03.rds")
    regions = c("Acinar","PanIN","PDAC")
    meta.df = Obj.list$meta %>% dplyr::filter(Group %in% regions)        
    plot.df = Obj.list$impute
    plot.df$pids = Obj.list$anno[row.names(plot.df),"pids"]
    plot.df$GeneNames = Obj.list$anno[row.names(plot.df),"GeneNames"]    
    
    dep.df = Obj.list$PDAC_stage %>% filter(log2FC>log2(2) & Pvalue<0.05)    
    plot.df = plot.df[dep.df$pid,]
    row.names(plot.df) = paste0(dep.df$pid,"_",dep.df$genename)
    plot_data_df = plot.df[,meta.df$SampleID]
    scale_data_df = t(apply(plot_data_df, 1, scale)) %>% as.data.frame()
    names(scale_data_df) =  names(plot_data_df)
    
    range_all <- range(c(-2, 2))    
    my_palette <- colorRampPalette(c("navy", "white","firebrick3"))(n=100)
    breaks = seq(range_all[1], range_all[2], length.out = 101)
  
    p1=pheatmap::pheatmap(scale_data_df,
                            scale = "none",show_rownames = F,show_colnames = T,
                            legend = T,border_color = NA,                        
                            color = my_palette, 
                            breaks = breaks,
                            cluster_cols = F,cluster_rows=T,
                            cellwidth = 25,cellheight = 0.3,fontsize = 14,fontsize_number = 10,
                            silent = F)
    pdf(file=paste0("write/Figure_4c.pdf"),width = 7,height = 8)
    print(p1)
    dev.off()
}

# Figure_4d, enrichment analysis of Sig. proteins for acinar, PanIN and PDAC regions
{
    rm(list=ls())
    Obj.list = readRDS(file = "sp_dataset_03.rds")
    dep.df = Obj.list$PDAC_stage %>% filter(Pvalue<0.05 & log2FC>log2(2))
    dep.df$ENTREZID = Obj.list$anno[dep.df$pid,"GeneID"]        

    formula_res <- compareCluster(data = dep.df,ENTREZID~region,fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                
    )
    formula_res_cutoff = formula_res
    formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue < 0.05,]
    formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)

    selected_GO = formula_res_cutoff@compareClusterResult
    selected_items = c("digestion","pancreatic juice secretion","intestinal cholesterol absorption","oxidative phosphorylation","lipid digestion",
                "Ras protein signal transduction","regulation of lymphocyte mediated immunity","regulation of stress fiber assembly","negative regulation of natural killer cell mediated cytotoxicity","regulation of actin filament bundle assembly",
                "positive regulation of response to wounding","barbed-end actin filament capping","positive regulation of wound healing","negative regulation of actin filament depolymerization","mitochondrial translation"
                )
    formula_res_cutoff@compareClusterResult = selected_GO %>% filter(Description %in% selected_items)
    
    p1=dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=14,color="log10P",size="count") + 
      theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
      scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30))
  
    pdf(file=paste0("write/Figure_4d.pdf"),width = 8,height = 6)
    print(p1)
    dev.off()
}

# for Figure_4g, heat map of the Sig. proteins after two-sided Student's t-test
{
    rm(list=ls())
    Obj.list = readRDS(file="sp_dataset_03.rds")
    regions = c("LN","IT")
    meta.df = Obj.list$meta %>% dplyr::filter(Group %in% regions)
        
    plot.df = Obj.list$impute[,meta.df$SampleID]
    dep.df = Obj.list$Lym_DEP %>% filter(abs(log2FC)>log2(2) & Pvalue<0.05)
            
    plot.df = plot.df[dep.df$pid,]
    row.names(plot.df) = paste0(dep.df$pid,"_",dep.df$genename)
    plot_data_df = plot.df
    scale_data_df = t(apply(plot_data_df, 1, scale)) %>% as.data.frame()
    names(scale_data_df) =  names(plot_data_df)
    
    range_all <- range(c(-2, 2))    
    my_palette <- colorRampPalette(c("navy", "white","firebrick3"))(n=100)
    breaks = seq(range_all[1], range_all[2], length.out = 101)
            
    p1=pheatmap::pheatmap(scale_data_df,
                            scale = "none",show_rownames = F,show_colnames = T,
                            legend = T,border_color = NA,                            
                            cluster_cols = F,cluster_rows=T,
                            color = my_palette,breaks=breaks,
                            cellwidth = 25,cellheight = 0.5,fontsize = 14,fontsize_number = 10,
                            silent = F)
    pdf(file=paste0("write/Figure_4g.pdf"),width = 7,height = 8)
    print(p1)
    dev.off()    
}

# for Figure_4h, enrichment analysis of Sig. proteins for IT and LN regions
{
    rm(list=ls())
    Obj.list = readRDS(file = "sp_dataset_03.rds")
    dep.df = Obj.list$Lym_DEP %>% filter(Pvalue<0.05 & abs(log2FC)>log2(2))
    dep.df$ENTREZID = Obj.list$anno[dep.df$pid,"GeneID"]
    dep.df$region = "LN"
    dep.df$region[which(dep.df$log2FC<0)] = "IT"

    formula_res <- compareCluster(data = dep.df,ENTREZID~region,fun="enrichGO", 
                                    OrgDb = org.Mm.eg.db,
                                    ont = "BP", pAdjustMethod = "BH",
                                    pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                    
    )
    formula_res_cutoff = formula_res    
    formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<0.05,] # save table

    # select cell-type specific enriched items for visualization
    selected_GO = formula_res_cutoff@compareClusterResult
    selected_items = c("T cell receptor signaling pathway","response to virus","regulation of lymphocyte mediated immunity","regulation of leukocyte cell-cell adhesion","B cell proliferation",        
        "response to extracellular stimulus","myeloid leukocyte differentiation","lymphocyte mediated immunity","extracellular matrix organization","cytokine production involved in immune response","antigen processing and presentation"
    )
        
    formula_res_cutoff@compareClusterResult = selected_GO %>% filter(Description %in% selected_items)    
    formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)
    p2=dotplot(formula_res_cutoff, label_format=50,showCategory=12,font.size=14,color="log10P",size="count") + 
        theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
        scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30))
    
    pdf(file=paste0("write/Fig4h_GOBP_Selected_compareCluster_",Sys.Date(),".pdf"),width = 8,height = 6)
    print(p2)
    dev.off()
}        
