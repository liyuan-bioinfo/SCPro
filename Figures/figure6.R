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
    Obj.list = readRDS(file = "ct_dataset_02.rds")
    temp.df = Obj.list$identified %>% t() %>% as.data.frame()
    
    meta.df = Obj.list$meta
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
    temp.df2 = Obj.list$quantified %>% t() %>% as.data.frame()
    
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
    plot.data["Total","protein_sum"]=dim(Obj.list$identified)[1]
    plot.data$CellType = row.names(plot.data)
    plot.data$protein_loca = "Identified"
    
    plot.data2["Total","protein_sum"]=dim(Obj.list$quantified)[1]
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
        axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1,color = c(levels(Obj.list$meta$group_color),"Black")),
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
    Obj.list = readRDS(file = "ct_dataset_02.rds")
    temp.df = Obj.list$min_10_impute_log2 %>% t() %>% as.data.frame()
    
    temp.df.anno = Obj.list$meta
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
