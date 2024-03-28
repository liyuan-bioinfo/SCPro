library(readxl)
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

setwd("")

# Figure_4c
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

# for Figure_4d
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

# for Figure_S9
{

    rm(list=ls())
    Obj_list = readRDS(file = "sp_dataset_03.rds")
    head(Obj_list$meta)

    meta_df = Obj_list$meta %>% filter(Group %in% c("Acinar","PanIN","PDAC"))
    meta_df$CellType = meta_df$Group        
    select_pid = c("Q7TPZ8","O09037","O88343","P07356")
    
    data_df = Obj_list$filter[select_pid,meta_df$SampleID]
    data_df$pid = row.names(data_df)        
    data_spread_df = data_df %>% tidyr::gather(key="SampleID",value="intensity",-pid)       
    data_spread_df$CellType = gsub(data_spread_df$SampleID, pattern="_.*$",replacement="")
    data_spread_df$CellType = factor(data_spread_df$CellType, levels = c("Acinar","PanIN","PDAC"))
    
    list_merge.list= list()    
    comparison_list = list(c("Acinar","PanIN"),c("Acinar","PDAC"),c("PanIN","PDAC"))    
    for (protein in select_pid){    
        temp_box.df = data_spread_df %>% dplyr::filter(pid == protein)
        title = paste0(select_pid,"_",Obj_list$anno[protein,"Gene"])
                
        p1=ggplot(data = temp_box.df,mapping=aes(x=CellType,y=intensity,fill=CellType),color="black")+
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
                axis.title = element_text(face="bold"),
                axis.text.y = element_text(colour="black", size = 9),
                axis.text.x = element_text(colour="black", size = 9),
                axis.line = element_line(size=0.5, colour = "black"))+
            labs(x="",y="log2 LFQ intensity",title = title)+
            ggpubr::stat_compare_means(comparisons = comparison_list, label="p.signif", paired = F,method="t.test")
        
        list_merge.list[[protein]] = p1
    }
    
    merge.plot = ggpubr::ggarrange(plotlist = list_merge.list,ncol = 4,nrow = 1)
    pdf(paste0("write/Figure_S9.pdf"),width = 12,height = 2.5,onefile = T)
    print(merge.plot)
    dev.off()    
}

# for Figure_4g
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

# for Figure_4h
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
    formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<0.05,] # save table

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
