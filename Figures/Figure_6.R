library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)

library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("")
# Figure_6c, heat map of Sig. PM proteins for 14 cell types
{
    rm(list=ls())
    Obj.list = readRDS(file = "ct_dataset_01.rds")
    pm.df = Obj.list$GO_mouse_pm %>% filter(loci=="pm")
    deg_df = Obj.list$dep_df[Obj.list$dep_df$pid %in% pm.df$pid,]    
    meta = Obj.list$meta
    
    plot.data = Obj.list$ND_impute_log2[unique(deg_df$pid),meta$SampleID]
    plot_data.median=data.frame(pid=row.names(plot.data))
    for(i in levels(meta$CellType)){
        temp_SampleID = meta[meta$CellType  %in% i,"SampleID"]
        plot_data.median[,i]=apply(plot.data[,temp_SampleID],1,median)
    }
    
    ## annotation
    annotation_col = data.frame(
        CellType = levels(meta$CellType)
    )
    row.names(annotation_col) = levels(meta$CellType)
    anno_colors = list(
        CellType = levels(meta$CellType_Color)
    )
    names(anno_colors$CellType) = levels(meta$CellType)
    
    row.names(plot_data.median) = plot_data.median$pid
    plot_data.median$pid = NULL
    
    range_all = range(c(-3, 3))    
    my_palette = colorRampPalette(c("navy", "white","firebrick3"))(n=100)
    breaks = seq(range_all[1], range_all[2], length.out = 101)
        
    scale_data_df = t(apply(plot_data.median, 1, scale)) %>% as.data.frame()
    names(scale_data_df) =  names(plot_data.median)

    # plot_data.median.t = plot_data.median %>% t() %>% as.data.frame()
    p1=pheatmap::pheatmap(scale_data_df,cluster_rows = F,show_rownames = F,cluster_cols = F,
            scale = "none",border_color = NA,legend = T,cellwidth = 20,cellheight = 3,                
            color = my_palette,breaks = breaks, 
            annotation_col = annotation_col,annotation_names_col = FALSE,
            show_colnames = F,legend_labels = F,annotation_legend = F,labels_col = NA,labels_row = NA,
            annotation_names_row  = FALSE,annotation_colors = anno_colors,silent=TRUE
    )

    pdf(file=paste0("write/Figure_6c_Pheatmap_Surface_",Sys.Date(),".pdf"),height = 10)
    print(p1)
    dev.off()


}

# Figure_6d, Line plot of top Sig. PM proteins for 14 cell types.
{

    
}



# Figure_6f, Volcano plot of Sig. proteins from two Treg subtypes
{
    rm(list=ls())
    obj_list = readRDS("ct_dataset_02.rds")
    dep_df = obj_list$dep_df
    logFC_cutoff <- log2(2)
    log10_P_Value_cutoff <- -log10(0.05)

    # add label
    known_marker = c("Cd8a","Cd8b","Cd4")
    dep_df$gene_plot = c()      
    temp_ind = dep_df$gene %in% known_marker
    dep_df$gene_plot[temp_ind] = dep_df[temp_ind,"gene"]
      
    # begin plot
    p1 = ggplot(data = dep_df,aes(x = logFC,y = `-log10_P.Value`))+
        geom_point(data = subset(dep_df,abs(logFC)<logFC_cutoff),size=2,
                   col = 'gray',alpha = 0.4)+
        geom_point(data = subset(dep_df,abs(`-log10_P.Value`)<log10_P_Value_cutoff & abs(logFC)>logFC_cutoff),
                   size=2,col = 'gray',alpha = 0.4)+
        geom_point(data = subset(dep_df,abs(`-log10_P.Value`)>log10_P_Value_cutoff & logFC>logFC_cutoff),
                   size=2,col = 'red',alpha = 1)+
        geom_point(data = subset(dep_df,abs(`-log10_P.Value`)>log10_P_Value_cutoff & logFC< -logFC_cutoff),
                   size=2,col = 'blue',alpha = 1)+
        geom_point(data = subset(dep_df,!is.na(gene_plot)),
                   size=2,colour="black",alpha = 0.4)+
        
        theme_bw()+
        labs(x='log2(fold change)',y='-log10(p-value)')+
        theme(legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              plot.title = element_text(hjust=0.5,size = 14),
              text = element_text(size=14),
              legend.position = 'none',
              axis.line = element_line(colour = "black"))+
        
        geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),lty = 3,col = 'black',lwd = 0.4)+
        geom_hline(yintercept = log10_P_Value_cutoff,lty = 3,col = 'black',lwd = 0.4)  +
        geom_text_repel(data = dep_df,min.segment.length=0, aes(label = gene_plot),size = 7,col = 'black') +
        ggtitle(paste("Volcano Plot of ",groupA," vs ",groupB,sep=""))
    
    # save plot
    pdf("Figure_6f_VolcanoPlot.pdf")
    print(p1)
    dev.off()    
}

# Figure_6g, Enrichment Analysis of Sig. Proteins from two Treg subtypes
{    
    rm(list=ls())
    obj_list = readRDS("ct_dataset_02.rds")
    dep_df = obj_list$dep_df    

    # set comparison group using Sig. proteins
    dep_df$group = ""
    dep_df$group[dep_df$P.Value<=0.05 & dep_df$logFC<=-1] = "Cd4+Cd25+Klrg1-"
    dep_df$group[dep_df$P.Value<=0.05 & dep_df$logFC>=1] = "Cd4+Cd25+Klrg1+"
    dep_df = dep_df %>% filter(group !="")

    # Uniprot-ID to Gene-ID
    sig.pid.tran <- bitr(row.names(dep_df), fromType = "UNIPROT", toType = c("ENSEMBL", "ENTREZID","SYMBOL"), 
                       OrgDb = org.Mm.eg.db)
    bg.pid.tran <- bitr(row.names(all.df), fromType = "UNIPROT", toType = c("ENSEMBL", "ENTREZID","SYMBOL"), 
                      OrgDb = org.Mm.eg.db)
    
    dep_df = dep_df[sig.pid.tran$UNIPROT,]
    dep_df$Entrez = sig.pid.tran$ENTREZID
    formula_res <- compareCluster(Entrez~group, data=dep_df, fun="enrichGO", 
                                OrgDb = org.Mm.eg.db,
                                ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 1,qvalueCutoff = 1,
                                universe = bg.pid.tran$ENTREZID
    )
    formula_res_cutoff = formula_res
    formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$pvalue<=0.05,]
    # save plot
    pdf(file="Fig6_g_Cd4+Cd25+Klrg1-vsCd4+Cd25+Klrg1+_Sig_Comparison_GO-BP_DotPlot.pdf",
      width=11,height=8)
    dotplot(formula_res,color="qvalue",showCategory = 10,label_format=50,font.size=16)
    dev.off()
    
}
