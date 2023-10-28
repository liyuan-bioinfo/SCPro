'---
title: "04_Fig5"
author: "Li-Yuan"
date: "2022-11-29"
output: github_document
---
'


## Protein Profiling
## Loading packages and setting global params
library(dplyr)
library(ggplot2) 
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(VennDiagram)

setwd("E:\\XYF（许燕芬）\\FACS_SISPROT\\workspace\\20221128_final\\04-surface_protein/")

## Fig5 B
## Protein Location Number across all Cell types
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221010.rds")
  temp.df = Obj.list$identified %>% t() %>% as.data.frame()
  all_mem.df = Obj.list$GO_mouse_pm
  meta.df = Obj.list$meta %>% select(SampleID,CellType) %>% unique()
  # plot.data = meta.df %>% select(CellType) %>% unique()
  # row.names(plot.data) = plot.data$CellType
  plot.data = data.frame()
  for(i in levels(Obj.list$meta$CellType)){
    temp_index = which(meta.df$CellType==i) # CellType index
    CellType.df = Obj.list$identified[,temp_index]
    CellType_identified.pid = row.names(CellType.df[which(rowSums(CellType.df)!=0),])
    temp_mem.df = all_mem.df[all_mem.df$pid %in% CellType_identified.pid,]
    sum.df=data.frame(CellType=i,membrane=table(temp_mem.df$loci)["membrane"], pm=table(temp_mem.df$loci)["pm"])
    plot.data = rbind(plot.data,sum.df)
  }
  
  plot.data["All",] = c("All",table(all_mem.df$loci)["membrane"],table(all_mem.df$loci)["pm"])
  plot_data.df = data.frame(CellType=c(plot.data$CellType,plot.data$CellType),
                            Number=c(plot.data$membrane,plot.data$pm),
                            Loci=c(rep("membrane",15),rep("pm",15)))
  plot_data.df$CellType = factor(plot_data.df$CellType, levels = c(levels(Obj.list$meta$CellType),"All"))
  plot_data.df$Loci = factor(plot_data.df$Loci, levels = c("membrane","pm"))
  plot_data.df$Number = as.numeric(plot_data.df$Number)
  plot_data.df$Number2 = plot_data.df$Number
  plot_data.df$Number2[which(plot_data.df$Loci=="pm")] = plot_data.df[which(plot_data.df$Loci=="pm"),"Number2"] * 3
  
  p1<-ggplot(data=plot_data.df, mapping=aes(x = CellType, y = Number,fill=Loci))+
    geom_bar(colour="black",stat="identity",position=position_dodge(0.9)) +
    theme_classic()+
    theme(
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          text=element_text(size = 14),
          axis.title = element_text(size=12),
          axis.text.y= element_text(size=13,colour = 'black'), # 设置y轴的数字大小
          axis.text.x = element_text(size=12,angle = 45, hjust = 1, vjust = 1,color = c(levels(Obj.list$meta$CellType_Color),"Black")),
          axis.ticks.x = element_blank(),axis.line.x.bottom = element_blank(),
          
          legend.position=c(0.5,0.8),legend.direction = "horizontal",
          legend.text = element_text(size=12),legend.background = element_blank()
    
          )+
    labs(x='',y='Identified membrane proteins',fill='')+
    theme(text=element_text(size=18,  color ='black')) + 
    scale_fill_grey(start = 0.5, end = 0.75) +
    scale_y_continuous(
                       expand = c(0,0),limits=c(0,1800),breaks=c(0,350,600,1000,1500,1800),label=c(0,350,600,1000,1500,1800))

    
    # scale_y_continuous(sec.axis = sec_axis(~. /3, name = "Identified plasma membrane proteins",breaks = c(0,200,400,600),labels = c(0,200,400,600)),
    #                    expand = c(0,0),limits=c(0,1800),breaks=c(0,500,1000,1500,1800),label=c(0,5000,1000,1500,1800))
    # 
  
    ##> save  
    pdf(p1,file=paste0("Fig5_B_Bar_Identified_Mem_genename_CellType_V2_",Sys.Date(),".pdf"),width = 5,height = 4)
    print(p1)
    dev.off()
    plot_data.df$Number2 = NULL
    write.csv(plot_data.df, file=paste0("Fig5_B_Bar_Identified_Mem_genename_CellType_",dim(plot_data.df)[1],"_",
                                    Sys.Date(),".csv"))
  
}

##DEG from one vs the rest, combine copynum,log2FC,protein function of surface
#'@return Obj.list$Fig5_deg
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  pm.df = Obj.list$GO_mouse_pm %>% filter(loci=="pm")#712
  deg_df = Obj.list$deg_Fig3[Obj.list$deg_Fig3$pid %in% pm.df$pid,]

  
  ##> copynum
  meta = Obj.list$meta
  deg_copynum.df = log2(Obj.list$identified_copynum[unique(deg_df$pid),meta$SampleID]+1)
  #mean or median
  deg_copynum.median=data.frame(pid=row.names(deg_copynum.df))
  for(i in levels(meta$CellType)){
    temp_SampleID = meta[meta$CellType  %in% i,"SampleID"]
    deg_copynum.median[,i]=apply(deg_copynum.df[,temp_SampleID],1,median)
  }
  
  deg_copynum_spread = deg_copynum.median %>% gather(key = "CellType",value = "median_copynum",-pid)
  
  temp.df = deg_df %>% merge(deg_copynum_spread,by="pid") %>% filter(CellType.x == CellType.y)
  
  temp.df = temp.df %>% merge(Obj.list$surface_class,by="pid",all.x=T,all.y=F)
  temp.df$CellType=temp.df$CellType.x
  temp.df = temp.df %>% select(pid,genename,CellType,Loca,logFC,median_copynum,class)
  
  temp.df$CellType = factor(temp.df$CellType,levels = levels(Obj.list$meta$CellType))
  temp.df = temp.df[order(temp.df$CellType),]
  write.csv(temp.df,file=paste0("Fi5_C_DEG_PM_",dim(temp.df)[1],"_",Sys.Date(),".csv"),row.names = F)

  Obj.list$deg_Fig5 = temp.df
  # saveRDS(Obj.list,file = "../FACS.MSFragger.Obj.20221010.rds")

}

## Fig5 C,pheatmap of deg pm proteins
{
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221128.rds")
  pm.df = Obj.list$GO_mouse_pm %>% filter(loci=="pm")
  deg_df = Obj.list$deg_Fig3[Obj.list$deg_Fig3$pid %in% pm.df$pid,]
  # cell_specific_pid = names(table(deg_df$pid))[-which(table(deg_df$pid)>1)]
  # deg_df = deg_df[deg_df$pid %in% cell_specific_pid,]
  write.csv(deg_df,file="surface.csv")
  
  meta = Obj.list$meta
  
  ##> plot
  plot.data = Obj.list$mann_imputate_log2[unique(deg_df$pid),meta$SampleID]
  #mean or median
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
  
  # plot_data.median.t = plot_data.median %>% t() %>% as.data.frame()
  pdf(file=paste0("Fig5_C_Pheatmap_Surface_deg_",Sys.Date(),".pdf"),height = 10)
  pheatmap(plot_data.median,cluster_rows = F,show_rownames = F,cluster_cols = F,
           scale = "row",border_color = NA,legend = F,cellwidth = 20,cellheight = 3,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(1000),
           annotation_col = annotation_col,annotation_names_col = FALSE,
           show_colnames = F,legend_labels = F,annotation_legend = F,labels_col = NA,labels_row = NA,
           annotation_names_row  = FALSE,annotation_colors = anno_colors
  )
  dev.off()
  plot_data.scale = apply(plot_data.median,1,scale) %>% t %>% as.data.frame()
  names(plot_data.scale) = names(plot_data.median)
  write.csv(plot_data.scale,file="pheatmap.csv")
  
}

## Fig5 S3
## Pie plot of all/deg pm proteins
{
  
  ##> pie plot for deg surface annotation category
  rm(list=ls())
  Obj.list = readRDS(file = "../FACS.MSFragger.Obj.20221010.rds")

  
  ##> pie plot for surface annotation category
  pm_anno.df = Obj.list$surface_class %>% group_by(class) %>% summarise(n=count(class))
  
  myPalette <- brewer.pal(8, "Set2") 
  
  # You can change the border of each area with the classical parameters:
  pdf(file=paste0("Fig5_S3_surface_anno_",Sys.Date(),".pdf"))
  pie(pm_anno.df$n$freq , labels = paste0(pm_anno.df$n$x,"(",pm_anno.df$n$freq,")"), border="white", col=myPalette,)
  dev.off()
  
  ##> pie plot for deg surface
  pm_anno.df = Obj.list$deg_Fig5 %>% group_by(class) %>% summarise(n=count(class))
  pdf(file=paste0("Fig5_S3_deg_surface_anno_",Sys.Date(),".pdf"))
  pie(pm_anno.df$n$freq , labels = paste0(pm_anno.df$n$x,"(",pm_anno.df$n$freq,")"), border="white", col=myPalette )
  dev.off()
    
  }

##Fig5 D-END
### Profile of select pm proteins with scores
{
  
  rm(list=ls())
  Obj.list = readRDS("../FACS.MSFragger.Obj.20221128.rds")
  
  df = Obj.list$min_10_imputate_log2
  rowmeans.df = data.frame(row.names = levels(Obj.list$meta$CellType))
  
  #row_means
  temp.mean.df = data.frame(row.names = row.names(Obj.list$filter))
  for(i in levels(Obj.list$meta$CellType)){
    temp_index = which(Obj.list$meta$CellType==i)  
    temp_df = Obj.list$min_10_imputate_log2[,temp_index]
    temp.mean.df[,i] = rowMeans(temp_df,na.rm = TRUE)
    # apply(plot.data[,temp_SampleID],1,mean)
  }
  #zscore
  temp.mean.zscore=apply(temp.mean.df,1,function(x){ 
    (x-mean(x))/sd(x)
  })
  temp.mean.zscore = temp.mean.zscore %>% t() %>% as.data.frame()
  Obj.list$min_10_imputate_log2_row_means_zscore = temp.mean.zscore
  run_comparision_line = function(temp_target_gene,temp_target_cell){
    # temp_target_cell = "PCC"
    # temp_target_gene = c("Epcam","Msln","Mmp7","Krt18")
    
    temp_bg_gene_meta = Obj.list$deg_Fig3 %>% dplyr::filter(CellType==temp_target_cell) %>% 
      dplyr::select(pid,genename) %>% unique()
    temp_bg_gene_meta$genename = Obj.list$annotation[temp_bg_gene_meta$pid,"genename"]
    
    temp_target_gene_meta = Obj.list$deg_Fig3 %>% dplyr::filter(CellType==temp_target_cell) %>% 
      dplyr::filter(genename %in% temp_target_gene )%>% 
      dplyr::select(pid,genename) %>% unique()
    row.names(temp_target_gene_meta) = temp_target_gene_meta$genename
    temp_target_gene_meta = temp_target_gene_meta[temp_target_gene,]
    
    temp.1.df = Obj.list$min_10_imputate_log2_row_means_zscore[temp_bg_gene_meta$pid,]
    temp.2.df = Obj.list$min_10_imputate_log2_row_means_zscore[temp_target_gene_meta$pid,]
    
    temp.plot.data.1 = temp.1.df %>% 
      gather(key = "CellType",value="zscore")
    temp.plot.data.1$genename = rep(temp_bg_gene_meta$genename,dim(temp.1.df)[2])
    
    temp.plot.data.2 = temp.2.df %>% 
      gather(key = "CellType",value="zscore")
    temp.plot.data.2$genename = rep(temp_target_gene_meta$genename,dim(temp.2.df)[2])
    temp.plot.data.2$genename = factor(temp.plot.data.2$genename,levels = temp_target_gene)
    
    temp.plot.color = c(brewer.pal(9, "YlOrRd")[9],brewer.pal(9, "YlGnBu")[9],
                        brewer.pal(9, "Greens")[9])#,brewer.pal(9, "BuPu")[9])
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
    temp.box.all = Obj.list$min_10_imputate_log2
    # temp.box.all = temp.box.all[c("Q99JW5","Q78IQ7","P55012","P09803"),]
    temp_target_meta = Obj.list$annotation[Obj.list$annotation$genename%in%temp_target_gene,]
    row.names(temp_target_meta) = temp_target_meta$genename
    temp_target_meta = temp_target_meta[temp_target_gene,]
    
    temp.box.all = temp.box.all[temp_target_meta$pid,]
    
    temp_color = c(brewer.pal(9, "YlOrRd")[9],brewer.pal(9, "YlGnBu")[9],
                   brewer.pal(9, "Greens")[9])#,brewer.pal(9, "BuPu")[9])
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
  temp_target_gene = c("Epcam","Cdh1","Itga6")
  temp_target_cell = "PCC"
  temp_plot.PCC = run_comparision_line(temp_target_gene, temp_target_cell)
  temp_plot_list.PCC = run_comparision_box(temp_target_gene)
  
  ##CAF
  temp_target_gene = c("Ptk7","Ror2","Pdgfrb")
  temp_target_cell = "CAF"
  temp_plot.CAF = run_comparision_line(temp_target_gene, temp_target_cell)
  temp_plot_list.CAF = run_comparision_box(temp_target_gene)
  
  ## Treg
  temp_target_gene = c("Tnfrsf18","Klrg1","Pdcd1")
  temp_target_cell = "Treg"
  temp_plot.Treg = run_comparision_line(temp_target_gene, temp_target_cell)
  
  temp_plot_list.Treg = run_comparision_box(temp_target_gene)
  
  ## DC
  temp_target_gene = c("Icam1","F11r","Itgae")
  temp_target_cell = "DC"
  temp_plot.DC = run_comparision_line(temp_target_gene, temp_target_cell)
  temp_plot_list.DC = run_comparision_box(temp_target_gene)
  
  p1=ggpubr::ggarrange(temp_plot.PCC,temp_plot_list.PCC[[1]],temp_plot_list.PCC[[2]],temp_plot_list.PCC[[3]],
                       temp_plot.CAF,temp_plot_list.CAF[[1]],temp_plot_list.CAF[[2]],temp_plot_list.CAF[[3]],
                       temp_plot.Treg,temp_plot_list.Treg[[1]],temp_plot_list.Treg[[2]],temp_plot_list.Treg[[3]],
                       temp_plot.DC,temp_plot_list.DC[[1]],temp_plot_list.DC[[2]],temp_plot_list.DC[[3]],
                    ncol = 4,nrow = 4
                    )
  pdf(file=paste0("Fig5_DEG_profile_",Sys.Date(),".pdf"),width = 10)
  print(p1)
  dev.off()
  
}


#------------
##NOT use
## Fig5 C, One vs the Rest
{
  run_OnevsRest = function(temp_meta, temp_df){
    
    # temp_meta = temp_submeta
    # temp_df = temp_subdf
    
    temp.output.list = list()
    for(i in levels(temp_meta$CellType)[-1]){
      temp_samples_A = temp_meta$SampleID[which(temp_meta$CellType==levels(temp_meta$CellType)[1])]
      temp_samples_B = temp_meta$SampleID[which(temp_meta$CellType==i)]
      temp_group_A = levels(temp_meta$CellType)[1]
      temp_group_B = i
      
      temp_meta.select = temp_meta %>% filter(CellType==temp_group_A | CellType==temp_group_B)
      temp_meta.select$CellType = factor(temp_meta.select$CellType, levels = c(temp_group_A,temp_group_B))
      temp.df.select = temp_df %>% select(temp_samples_A,temp_samples_B) #df
      
      # temp.design <- model.matrix(~0+factor(temp_meta$CellType))
      temp.design <- model.matrix(~temp_meta.select$CellType)
      colnames(temp.design) <- levels(temp_meta.select$CellType)
      rownames(temp.design) <- temp_meta.select$SampleID
      
      
      fit <- lmFit(temp.df.select, temp.design)
      fit <- eBayes(fit, trend=TRUE)
      temp.result.limma <- topTable(fit,  coef=2,n=Inf)
      
      temp.df.select = temp.df.select[row.names(temp.result.limma),]
      temp.df.select$P.Value = temp.result.limma$P.Value
      temp.df.select$fdr = temp.result.limma$adj.P.Val
      temp.df.select$logFC = temp.result.limma$logFC
      
      # print(temp.df.select[c("P25918"),])
      
      temp.df.select.sig = temp.df.select %>% filter(P.Value<0.05 & logFC < 0)
      # temp.df.select.sig = temp.df.select.sig[order(abs(temp.df.select.sig$logFC),decreasing = TRUE),]
      
      temp.output.list[[paste(temp_group_A,"-",temp_group_B,sep="")]] = temp.df.select.sig
    }
    ## all-diff gene
    temp.x=row.names(temp.output.list[[1]])
    if(length(temp.x)==0){
      return(0)
    }
    
    # if(length(temp.output.list)==1){
    #   # temp.x = intersect(temp.x,row.names(temp.output.list[[1]]))
    # }else{
    for(i in (1:length(temp.output.list))){
      temp.x = intersect(temp.x,row.names(temp.output.list[[i]]))
      print(length(temp.x))
    }
    # }
    if(length(temp.x)==0){
      return(NULL)
    }
    
    
    ## average FC
    temp.df.sig = data.frame(pid=temp.x,stringsAsFactors = F)
    temp.logFC = temp.output.list[[1]]
    temp.logFC = temp.logFC[temp.df.sig$pid,"logFC"]
    
    if(length(temp.output.list)==1){
      
    }else{
      
      for(i in (2:length(temp.output.list))){
        t = temp.output.list[[i]][temp.df.sig$pid,"logFC"]
        temp.logFC = temp.logFC+t
        # print(length(x))
      }
    }
    temp.df.sig$logFC = abs(temp.logFC)
    temp.df.sig = temp.df.sig[order(temp.df.sig$logFC,decreasing = T),]
    temp.df.sig$CellType=levels(temp_meta$CellType)[1]
    return(temp.df.sig)
  }
  ### PCC vs Rest
  meta <- Obj.list$meta
  
  temp_group = meta %>% dplyr::select(CellType) %>% unique() 
  row.names(temp_group) = temp_group$CellType
  temp_group$CellType = as.character(temp_group$CellType)
  
  temp_OnevsRest_list = list()
  
  ##> PCC
  temp_subgroup = unique(c("PCC",temp_group$CellType))
  temp_OnevsRest_list[["PCC"]] = factor(temp_subgroup, levels = temp_subgroup)
  
  ##>CAF
  temp_subgroup = unique(c("apCAF",temp_group$CellType))
  temp_OnevsRest_list[["apCAF"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("myCAF",temp_group$CellType))
  temp_OnevsRest_list[["myCAF"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("iCAF",temp_group$CellType))
  temp_OnevsRest_list[["iCAF"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("CAF",temp_group$CellType))
  temp_OnevsRest_list[["CAF"]] = factor(temp_subgroup, levels = temp_subgroup)
  
  ##>Lymphoid
  temp_subgroup = unique(c("T4",temp_group$CellType))
  temp_OnevsRest_list[["T4"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("T8",temp_group$CellType))
  temp_OnevsRest_list[["T8"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("Treg",temp_group$CellType))
  temp_OnevsRest_list[["Treg"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("B",temp_group$CellType))
  temp_OnevsRest_list[["B"]] = factor(temp_subgroup, levels = temp_subgroup)
  
  ##> MYE
  temp_subgroup = unique(c("MYE",temp_group$CellType))
  temp_OnevsRest_list[["MYE"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("NEU",temp_group$CellType))
  temp_OnevsRest_list[["NEU"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("MO",temp_group$CellType))
  temp_OnevsRest_list[["MO"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("MAC",temp_group$CellType))
  temp_OnevsRest_list[["MAC"]] = factor(temp_subgroup, levels = temp_subgroup)
  temp_subgroup = unique(c("DC",temp_group$CellType))
  temp_OnevsRest_list[["DC"]] = factor(temp_subgroup, levels = temp_subgroup)    
  
  temp_out.df = data.frame()
  for(i in temp_OnevsRest_list){
    temp_submeta = meta %>% filter(CellType %in% i)
    temp_submeta$CellType = factor(x=temp_submeta$CellType, levels = i)
    temp_submeta = temp_submeta[order(temp_submeta$CellType),]
    
    temp_subdf = Obj.list$min_10_imputate_log2 %>% 
      select(temp_submeta$SampleID)
    temp.out = run_OnevsRest(temp_submeta, temp_subdf)
    if(!is.null(temp.out)){
      temp_out.df = rbind(temp_out.df,temp.out)  
    }
    
    print(i)
    # break
  }
  
  temp_out.df$genename = Obj.list$annotation[temp_out.df$pid,"genename"]
  pm.df = Obj.list$GO_mouse_pm %>% filter(loci=="pm")
  temp_out.df[temp_out.df$pid %in% pm.df$pid,"Loca"]="pm"
  write.csv(temp_out.df,file=paste0("Fig5_OnevsRest_DEG_Min10Imputate_V2_",dim(temp_out.df)[1],"_",Sys.Date(),".csv"),row.names = F,na = "")
  
  temp_out.df$CellType = factor(temp_out.df$CellType,levels=levels(Obj.list$meta$CellType))
  temp_out.df = temp_out.df[order(temp_out.df$CellType),]
  
  # Obj.list$deg_Fig5 = temp_out.df
  # saveRDS(Obj.list,file = "../FACS.MSFragger.Obj.20221010.rds")
}
