library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)

library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("")


# Bar plot showing the protein numbers of quantified and identified proteins 
{
    rm(list=ls())
    obj_list = readRDS(file = "ct_dataset_02.rds")
    temp.df = obj_list$identified %>% t() %>% as.data.frame()
    
    meta.df = obj_list$meta
    meta.df$SampleID = meta.df$sample_id
    meta.df$CellType = meta.df$group
    meta.df = meta.df %>% select(SampleID,CellType) %>% unique()
    # row.names(plot.data) = plot.data$group
    plot.data = meta.df %>% select(CellType) %>% unique()
    row.names(plot.data) = plot.data$CellType
    
    temp_sum = c()
    for(i in levels(plot.data$CellType)){
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
    plot.data2 = meta.df %>% select(CellType) %>% unique()
    row.names(plot.data2) = plot.data2$CellType
    temp.df2 = obj_list$quantified %>% t() %>% as.data.frame()
    
    temp_sum = c()
    for(i in levels(plot.data2$CellType)){
      temp_index = which(meta.df$CellType==i)  
      # t = unique(which(temp.df[temp_index,]!=0))
      tt = temp.df2[temp_index,]
      ttt = apply(tt,2,function(x){
        return(sum(x,na.rm = TRUE))
      })
      temp_len = length(which(ttt >0))
      temp_sum = c(temp_sum,temp_len)
      
    }
    plot.data2$protein_sum =temp_sum
    
    #combine
    plot.data["Total","protein_sum"]=dim(obj_list$identified)[1]
    plot.data$CellType = row.names(plot.data)
    plot.data$protein_loca = "Identified"
    
    plot.data2["Total","protein_sum"]=dim(obj_list$quantified)[1]
    plot.data2$CellType = row.names(plot.data2)
    plot.data2$protein_loca = "Quantified"
    
    plot.data3 = plot.data
    plot.data3$protein_sum = plot.data$protein_sum - plot.data2$protein_sum
    plot.data.bar = rbind(plot.data3,plot.data2)
    plot.data.bar$CellType = factor(x = plot.data.bar$CellType,levels = c(levels(meta.df$CellType),"Total"))
    
    temp.p1 <- ggplot(plot.data.bar, aes(x=CellType,y=protein_sum,fill=protein_loca)) + 
      geom_bar(stat = "identity",width=0.8)+
      scale_fill_grey(start = 0.75, end = 0.5) +
      theme_classic()+
      theme(        
        plot.title = element_text(size=16,hjust = 0.5),
        axis.title.y.left  = element_text(size=16),
        axis.text = element_text(size=12),
        axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1,color = c(levels(obj_list$meta$group_color),"Black")),
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
              file = paste0("Fig2_C_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".csv"),
              fileEncoding = "UTF-8",row.names = F,quote = TRUE)
    
    pdf(paste0("Fig2_C_ProteinNumber_PG_IdentifiedvsQuantified_Barplot_",Sys.Date(),".pdf"),
        width=6,height = 4)
    print(temp.p1)
    dev.off()
  }   

# PCA plot showing the 8 sub-type of Treg
{
    rm(list=ls())
    obj_list = readRDS(file = "ct_dataset_02.rds")
    temp.df = obj_list$min_10_impute_log2 %>% t() %>% as.data.frame()
    
    temp.df.anno = obj_list$meta
    temp.df.anno = temp.df.anno[row.names(temp.df),]
    temp.df.anno$sample_id = c(
      paste0("CD4+CD25+Klrg1-_",1:3),paste0("CD4+CD25+Klrg1+_",1:3),
      paste0("CD4+CD25-_",1:3),paste0("CD4+CD25+_",1:3),
      paste0("CD4+Klrg1-_",1:3),paste0("CD4+Klrg1+_",1:3),
      paste0("CD8+Klrg1-_",1:3),paste0("CD8+Klrg1+_",1:3)
    )
    temp.df.anno$group =  gsub(temp.df.anno$sample_id,replacement = "",pattern = "_.*$")
    temp.df.anno$group = factor(temp.df.anno$group,levels = (
      c("CD4+CD25+Klrg1-","CD4+CD25+Klrg1+",
      "CD4+CD25-","CD4+CD25+",
      "CD4+Klrg1-","CD4+Klrg1+",
      "CD8+Klrg1-","CD8+Klrg1+")
      
    ))
    temp.df.anno$shape = factor(      
      c(rep("CD4+CD25+Klrg1",6),
        rep("CD4+CD25",6),
        rep("CD4+Klrg1",6),
        rep("CD8+Klrg1",6)),
      
      levels = (
      c("CD4+CD25+Klrg1",
        "CD4+CD25",
        "CD4+Klrg1",
        "CD8+Klrg1")
      
    ))
    pca_result <- prcomp(temp.df,scale=T,center = T)
    
    df1 = as.data.frame(pca_result$x)
    summ1 <- summary(pca_result)
    xlab1 <- paste0("PC1 (",round(summ1$importance[2,1]*100,2),"%)")
    ylab1 <- paste0("PC2 (",round(summ1$importance[2,2]*100,2),"%)")

    # begion plot
    p1=ggplot(data = df1,aes(x = PC1,y = PC2,color = temp.df.anno$group,
                             shape=temp.df.anno$shape
                             ))+

      geom_point(size = 3.5)+
      # geom_text_repel(data=df1,aes(x=PC1,y=PC2,color=temp.df.anno$group,label=temp.df.anno$sample_id))+
      labs(x = xlab1,y = ylab1,color = "",title = "",shape="")+
      guides(fill = "none")+
      theme_bw()+
      # scale_fill_manual(values = levels(temp.df.anno$))+
      scale_colour_manual(values = levels(temp.df.anno$group_color))+
      theme(plot.background = element_blank(),legend.background = element_blank(),
            panel.background = element_blank(),panel.grid = element_blank(),
            axis.text = element_text(size = 12),axis.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
    
    # save plot
    pdf(file=paste("1_ProteinExpression_Cd25vsKlrg1_All/Fig6_C_PCA_Minimputed.",Sys.Date(),".pdf",sep=""),
        width = 8,height=5)
    print(p1)
    dev.off()
    
}

# Volcano plot
{
    rm(list=ls())
    obj_list = readRDS("Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")
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
    pdf("VolcanoPlot.pdf")
    print(p1)
    dev.off()    
}

# Enrichment Analysis of Sig. Proteins
{    
    rm(list=ls())
    obj_list = readRDS("Validation_Cd25_Klrg1_All_Obj_2022-11-22.rds")
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
    
    # save plot
    pdf(file="Fig6_F_Cd4+Cd25+Klrg1-vsCd4+Cd25+Klrg1+_Sig_Comparison_GO-BP_DotPlot_2022-07-28.pdf",
      width=11,height=8)
    dotplot(formula_res,color="pvalue",showCategory = 10,label_format=50,font.size=16)
    dev.off()
    
}
