library(readxl)
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggplot2)

setwd("")

# Figure_4c
{
    rm(list=ls())
    Obj.list = readRDS(file="sp_dataset_02.rdata")
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

# for Figure-4d
{
    rm(list=ls())

    Obj.list = readRDS(file = "write/sp_dataset_02.rdata")
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

    
    p1=dotplot(formula_res_cutoff, label_format=50,showCategory=5,font.size=14,color="log10P",size="count") + 
      theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
      scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")[3:6])(30))
  
    pdf(file=paste0("write/Figure_4d.pdf"),width = 8,height = 6)
    print(p1)
    dev.off()
}

    ## Box plto
    # for Fig_S9
    {
    
        rm(list=ls())
        Obj_list = readRDS(file = "sp_dataset_02.rdata")
        names(Obj_list)
        head(Obj_list$meta)

        meta_df = Obj_list$meta %>% filter(Group %in% c("Acinar","PanIN","PDAC"))# 12 * 4
        meta_df$CellType = meta_df$Group        

        select_pid = c("Q7TPZ8","O09037","O88343","P07356")
        
        data_df = Obj_list$filter[select_pid,meta_df$SampleID] # use impute data        
        data_df$pid = row.names(data_df)        
        data_spread_df = data_df %>% tidyr::gather(key="SampleID",value="intensity",-pid)       
        data_spread_df$CellType = gsub(data_spread_df$SampleID, pattern="_.*$",replacement="")
        data_spread_df$CellType = factor(data_spread_df$CellType, levels = c("Acinar","PanIN","PDAC"))
        
                        
        
        list_merge.list= list()
        
        comparison_list = list(c("Acinar","PanIN"),c("Acinar","PDAC"),c("PanIN","PDAC"))
        
        for (protein in select_pid){
        # # print(protein)
        
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
                labs(x="",y="log2 LFQ intensity",title = title)+
                ggpubr::stat_compare_means(comparisons = comparison_list, label="p.signif", paired = F,method="t.test")
            
            list_merge.list[[protein]] = p1
        }
        
        merge.plot = ggpubr::ggarrange(plotlist = list_merge.list,ncol = 4,nrow = 1)
        pdf(paste0("write/Fig_S9_PDAC_Sig_Quantified_Box_",Sys.Date(),".pdf"),width = 12,height = 2.5,onefile = T)
        print(merge.plot)
        dev.off()
        
        
        ## save data
        write.csv(data_df,
                file = paste0("write/Sig_S9_PDAC_Sig_Quantified_Box_",Sys.Date(),".csv"),row.names = T)
    
    
    }
    

}

# IT and LN
{

    ## T.test
    {
        rm(list=ls())
        Obj.list = readRDS(file="write/sp_dataset_02.rdata")
        regions = c("LN","IT")
        meta.df = Obj.list$meta %>% dplyr::filter(Group %in% regions)
            
        meta.df$Regions = factor(meta.df$Group,levels = regions)
        meta.df = meta.df %>% arrange(Regions)
        input.df = Obj.list$filter[,meta.df$SampleID]
        input.df = input.df[-which(rowSums(input.df)==0),]
        
        Pvalue = c()
        log2FC = c()
        for(i in 1:dim(input.df)[1]){
            pid.df = input.df[i,] %>% t() %>% as.data.frame()
            names(pid.df) = "pid"
            pid.df$SampleId = row.names(pid.df)
            pid.df$Region = meta.df$Regions
            p.value = t.test(data=pid.df,pid~Region)$p.value
            Pvalue = c(Pvalue,p.value)
            
            temp_FC = median(pid.df$pid[1:3]) - median(pid.df$pid[4:6])
            log2FC = c(log2FC, temp_FC)
            
        }
        input.df$Pvalue = Pvalue
        input.df$fdr = p.adjust(Pvalue,method = "BH")
        input.df$log2FC = log2FC
        
        input.df$pid = Obj.list$anno[row.names(input.df),"pid"]
        input.df$genename = Obj.list$anno[row.names(input.df),"Gene"]
        
        # Obj.list$Lym_DEP = input.df
        # saveRDS(Obj.list,file = "write/sp_dataset_02.rdata")
        write.csv(input.df,file="write/Fig4g_ttest_FC2_table1.csv",row.names = F)
    
    
    }

    # Pheatmap
    {
        rm(list=ls())
        Obj.list = readRDS(file="sp_dataset_02.rdata")
        regions = c("LN","IT")
        meta.df = Obj.list$meta %>% dplyr::filter(Group %in% regions)
            
        plot.df = Obj.list$impute[,meta.df$SampleID]
        dep.df = Obj.list$Lym_DEP %>% filter(abs(log2FC)>log2(2) & Pvalue<0.05)
                
        plot.df = plot.df[dep.df$pid,]
        row.names(plot.df) = paste0(dep.df$pid,"_",dep.df$genename)
        
        range_all <- range(c(-2, 2))    
        my_palette <- colorRampPalette(c("navy", "white","firebrick3"))(n=100)
        breaks = seq(range_all[1], range_all[2], length.out = 101)
            
        plot_data_df = plot.df
        scale_data_df = t(apply(plot_data_df, 1, scale)) %>% as.data.frame()
        names(scale_data_df) =  names(plot_data_df)
        
        p1=pheatmap::pheatmap(scale_data_df,
                                scale = "none",show_rownames = F,show_colnames = T,
                                legend = T,border_color = NA,
                                
                                cluster_cols = F,cluster_rows=T,
                                color = my_palette,breaks=breaks,
                                cellwidth = 25,cellheight = 0.5,fontsize = 14,fontsize_number = 10,
                                silent = F)
        pdf(file=paste0("write/Fig4g_dep_proteins_pheatmap_2FC_",Sys.Date(),".pdf"),width = 7,height = 8)
        print(p1)
        dev.off()    
    }

    {
    ## GO-Anno
    {
        rm(list=ls())

        Obj.list = readRDS(file = "write/sp_dataset_02.rdata")
        dep.df = Obj.list$Lym_DEP %>% filter(Pvalue<0.05 & abs(log2FC)>log2(2))
        dep.df$ENTREZID = Obj.list$anno[dep.df$pid,"GeneID"]
        dep.df$region = "LN"
        dep.df$region[which(dep.df$log2FC<0)] = "IT"
        # tran.df = bitr(dep.df$pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        # dep.df = dep.df %>% merge(tran.df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)

        # bg_pid = row.names(Obj.list$impute)
        # bg_pid.df = bitr(bg_pid, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

        formula_res <- compareCluster(data = dep.df,ENTREZID~region,fun="enrichGO", 
                                        OrgDb = org.Mm.eg.db,
                                        ont = "BP", pAdjustMethod = "BH",
                                        pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                        
        )
        formula_res_cutoff = formula_res
        
        formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<0.05,]
        # write.csv(x = formula_res_cutoff@compareClusterResult, file=paste0("write/04_LN_IT/pvalue_sig_GOBP_",Sys.Date(),".csv"))  

        # Load 2023-05-07 GO-BP
        selected_GO = read.csv("/aaa/zihanwu/yyyli2/projectx_xu_NC2023/datasets/dataset_03/write/04_LN_IT/pvalue_sig_GOBP_2023-05-07.csv", row.names=1)
        # selected_GO = formula_res_cutoff@compareClusterResult
        # selected_GO = formula_res_cutoff@compareClusterResult
        selected_items = c("T cell receptor signaling pathway","response to virus","regulation of lymphocyte mediated immunity","regulation of leukocyte cell-cell adhesion","B cell proliferation",        
        "response to extracellular stimulus","myeloid leukocyte differentiation","lymphocyte mediated immunity","extracellular matrix organization","cytokine production involved in immune response","antigen processing and presentation"

        )
        
        selected_GO = selected_GO %>% filter(Description %in% selected_items)
        formula_res_cutoff@compareClusterResult = selected_GO %>% filter(p.adjust < 0.05)
        
        formula_res_cutoff@compareClusterResult$log10P = -log10(formula_res_cutoff@compareClusterResult$pvalue)
        p2=dotplot(formula_res_cutoff, label_format=50,showCategory=12,font.size=14,color="log10P",size="count") + theme(panel.grid = element_blank(),axis.ticks.y = element_blank()) +
        scale_colour_gradientn(colours=colorRampPalette(RColorBrewer::brewer.pal(n = 7, 
                                                                                name = "YlOrRd")[3:6])(30))
        pdf(file=paste0("write/Fig4h_GOBP_Selected_compareCluster_",Sys.Date(),".pdf"),width = 8,height = 6)
        print(p2)
        dev.off()
    }        
    }
}


# Fig Legends
{
    rm(list=ls())
    # Fig lengends, -2 to 2           
    my_palette <- c(
        c(colorRampPalette(c("navy", "white"))(50)),
        c(colorRampPalette(c("white", "firebrick3"))(50))
    )        
    df <- data.frame(matrix(nrow = 10, ncol = 2))
    df[] <- rnorm(20)    
    p1=ggplot(data = df, aes(x = X1, y = X2, fill = X2)) + geom_tile() + theme_void()+#
        scale_fill_gradientn(colors=my_palette, limits = c(-2, 2), breaks = c(-2,0,2),guide = guide_colorbar(barwidth = 5, barheight = 1, 
                                                                                        ticks = TRUE, ticks.colour="black",ticks.linewidth=1/.pt,direction="horizontal",#"vertical",                                                                                            
                                                                                        label.theme = element_text(size = 8),
                                                                                        draw.ulim = TRUE, 
                                                                                        draw.llim = TRUE, 
                                                                                        draw.separator = TRUE,
                                                                                        
                                                                                        title.position = "top"
                                                                                        )        
        )

    p2=ggplot(data = df, aes(x = X1, y = X2, fill = X2)) + geom_tile() + theme_void()+#
        scale_fill_gradientn(colors=my_palette, limits = c(-2, 2), breaks = c(-2,0,2),guide = guide_colorbar(barwidth = 1, barheight = 5, 
                                                                                        ticks = TRUE, ticks.colour="black",ticks.linewidth=1/.pt,direction="vertical",#"vertical",                                                                                            
                                                                                        label.theme = element_text(size = 8),
                                                                                        draw.ulim = TRUE, 
                                                                                        draw.llim = TRUE, 
                                                                                        draw.separator = TRUE,
                                                                                        
                                                                                        title.position = "top"
                                                                                        )        
        )        

    pdf(file=paste0("write/Fig4_Fig_Legends_zscore_Pheatmap_",Sys.Date(),".pdf"))
    print(p1)
    print(p2)
    dev.off()  

}
