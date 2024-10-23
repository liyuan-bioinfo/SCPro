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
    rm(list=ls())
    Obj.list = readRDS("ct_dataset_01.rds")
    
    df = Obj.list$Mini_impute_log2
    rowmeans.df = data.frame(row.names = levels(Obj.list$meta$CellType))

    temp.mean.df = data.frame(row.names = row.names(Obj.list$filter))
    for(i in levels(Obj.list$meta$CellType)){
        temp_index = which(Obj.list$meta$CellType==i)  
        temp_df = Obj.list$Mini_impute_log2[,temp_index]
        temp.mean.df[,i] = rowMeans(temp_df,na.rm = TRUE)
        # apply(plot.data[,temp_SampleID],1,mean)
    }
    #zscore
    temp.mean.zscore=apply(temp.mean.df,1,function(x){ 
        (x-mean(x))/sd(x)
    })
    temp.mean.zscore = temp.mean.zscore %>% t() %>% as.data.frame()
    Obj.list$Mini_impute_scale = temp.mean.zscore
    run_comparision_line = function(temp_target_gene,temp_target_cell){
        # temp_target_cell = "PCC"
        # temp_target_gene = c("Epcam","Msln","Mmp7","Krt18")
        
        temp_bg_gene_meta = Obj.list$dep_df %>% dplyr::filter(CellType==temp_target_cell) %>% 
        dplyr::select(pid,genename) %>% unique()
        temp_bg_gene_meta$genename = Obj.list$annotation[temp_bg_gene_meta$pid,"genename"]
        
        temp_target_gene_meta = Obj.list$dep_df %>% dplyr::filter(CellType==temp_target_cell) %>% 
        dplyr::filter(genename %in% temp_target_gene )%>% 
        dplyr::select(pid,genename) %>% unique()
        row.names(temp_target_gene_meta) = temp_target_gene_meta$genename
        temp_target_gene_meta = temp_target_gene_meta[temp_target_gene,]
        
        temp.1.df = Obj.list$Mini_impute_scale[temp_bg_gene_meta$pid,]
        temp.2.df = Obj.list$Mini_impute_scale[temp_target_gene_meta$pid,]
        
        temp.plot.data.1 = temp.1.df %>% 
        gather(key = "CellType",value="zscore")
        temp.plot.data.1$genename = rep(temp_bg_gene_meta$genename,dim(temp.1.df)[2])
        
        temp.plot.data.2 = temp.2.df %>% 
        gather(key = "CellType",value="zscore")
        temp.plot.data.2$genename = rep(temp_target_gene_meta$genename,dim(temp.2.df)[2])
        temp.plot.data.2$genename = factor(temp.plot.data.2$genename,levels = temp_target_gene)
        
        temp.plot.color = c(brewer.pal(9, "YlOrRd")[9],brewer.pal(9, "YlGnBu")[9])#,
                            # brewer.pal(9, "Greens")[9])#,brewer.pal(9, "BuPu")[9])
        temp.plot.data.2$color = rep(temp.plot.color, dim(temp.2.df)[2])
        temp.plot.data.2$color = factor(temp.plot.data.2$color,levels = temp.plot.color)
        
        temp.plot.data.1$CellType = factor(temp.plot.data.1$CellType,levels=levels(Obj.list$meta$CellType))
        temp.plot.data.2$CellType = factor(temp.plot.data.2$CellType,levels=levels(Obj.list$meta$CellType))
        
        temp.plot = 
        ggplot() +
        geom_line(data=temp.plot.data.1,aes(x=CellType, y=zscore, group=genename),size=1,
                    color=c("darkgrey"))+
        geom_line(data=temp.plot.data.2,aes(x=CellType, y=zscore, group=genename, colour=genename),size=1.5,show.legend = T  
                    
        ) +
        # geom_text_repel( data=temp.line.anno.data, aes(x=rank,y=value,color=name, label = marker), hjust = -3, vjust = -1, size = 5.0,show.legend = F,min.segment.length = 1 ) +
        geom_point(data=temp.plot.data.2,aes(x=CellType,y=zscore,color=color),color=temp.plot.data.2$color,size=2,show.legend = F)+
        xlab("")+ylab("Z-score")+
        scale_color_manual(values =levels(temp.plot.data.2$color))+
        theme(axis.title = element_text(size = 12),axis.title.y = element_text(size = 12),legend.position = "right")+
        theme_classic()+
        theme(
            # panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background = element_blank(),
            legend.key=element_blank(),
            # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=10),
            axis.text.x = element_text(size=10,angle = 45, hjust = 1, vjust = 1,color = c(levels(Obj.list$meta$CellType_Color))),
            legend.position=c(0.5,0.85),legend.title = element_text(size = 5),
            legend.text = element_text(size = 8),
            legend.direction = "horizontal",legend.background = element_blank()
        )+
        guides(colour=guide_legend(title=""))
        return(temp.plot)
    }
    run_comparision_box = function(temp_target_gene){
        temp.box.all = Obj.list$Mini_impute_log2
        # temp.box.all = temp.box.all[c("Q99JW5","Q78IQ7","P55012","P09803"),]
        temp_target_meta = Obj.list$annotation[Obj.list$annotation$genename%in%temp_target_gene,]
        row.names(temp_target_meta) = temp_target_meta$genename
        temp_target_meta = temp_target_meta[temp_target_gene,]
        
        temp.box.all = temp.box.all[temp_target_meta$pid,]
        
        temp_color = c(brewer.pal(9, "YlOrRd")[9],brewer.pal(9, "YlGnBu")[9])#,
                    # brewer.pal(9, "Greens")[9])#,brewer.pal(9, "BuPu")[9])
        temp_plot_list = list()
        for(i in (1:dim(temp.box.all)[1])){
        
        temp.pid = row.names(temp.box.all[i,])
        temp.genename = Obj.list$annotation[temp.pid,"genename"]
        temp.df.gather = temp.box.all[i,] %>% gather(key="cell_type",value="value")
        
        temp.df.gather$name = gsub(temp.df.gather$cell_type,replacement = "",pattern = "_.*$")
        
        temp_plot_list[[i]] = ggpubr::ggboxplot(temp.df.gather, x = "name", y = "value",color=temp_color[i],
                                                # palette = "jco",
                                                add="jitter", add.params=list(size=1,fill=temp_color[i]),bxp.errorbar = TRUE,
                                                # order = c("PCC","CAF","myCAF","apCAF","iCAF",
                                                #           "T4","T8","Treg","B","MYE","DC","MAC",
                                                #           "MONO","NEU")
                                                order=levels(Obj.list$meta$CellType)
        )+ 
            # ggpubr::stat_compare_means(comparisons = my_comparisons,method = "t.test", 
            #                            method.args = list(var.equal = TRUE), 
            #                            label = "p.signif",hide.ns = TRUE,vjust = 0.5)+
            xlab("")+
            ylab("log2 LFQ intensity")+
            theme(legend.position = "none")+ggtitle(paste(temp.pid,temp.genename,sep="_"))+
            theme(plot.title = element_text(hjust = 0.5),
                text = element_text(size=10),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                axis.title.x=element_blank()
            )
        
        }
        return(temp_plot_list)
    }
    
        
    ## PCC
    temp_target_gene = c("Epcam","Cdh1")
    temp_target_cell = "PCC"
    temp_plot.PCC = run_comparision_line(temp_target_gene, temp_target_cell)
    temp_plot_list.PCC = run_comparision_box(temp_target_gene)
    
    ##CAF
    temp_target_gene = c("Ptk7","Ror2")
    temp_target_cell = "CAF"
    temp_plot.CAF = run_comparision_line(temp_target_gene, temp_target_cell)
    temp_plot_list.CAF = run_comparision_box(temp_target_gene)
    
    ## Treg
    temp_target_gene = c("Tnfrsf18","Klrg1")
    temp_target_cell = "Treg"
    temp_plot.Treg = run_comparision_line(temp_target_gene, temp_target_cell)
    
    temp_plot_list.Treg = run_comparision_box(temp_target_gene)
    
    ## DC
    temp_target_gene = c("Icam1","F11r")
    temp_target_cell = "DC"
    temp_plot.DC = run_comparision_line(temp_target_gene, temp_target_cell)
    temp_plot_list.DC = run_comparision_box(temp_target_gene)
    
    p1=ggpubr::ggarrange(temp_plot.PCC,temp_plot_list.PCC[[1]],temp_plot_list.PCC[[2]],
                        temp_plot.CAF,temp_plot_list.CAF[[1]],temp_plot_list.CAF[[2]],
                        temp_plot.Treg,temp_plot_list.Treg[[1]],temp_plot_list.Treg[[2]],
                        temp_plot.DC,temp_plot_list.DC[[1]],temp_plot_list.DC[[2]],
                        ncol = 3,nrow = 4
                        )
    pdf(file=paste0("write/Figure_6d_DEP_profile_",Sys.Date(),".pdf"),width = 10)
    print(p1)
    dev.off()

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
