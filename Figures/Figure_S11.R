library(readxl)
library(dplyr)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

setwd("")

# for Figure_S11 a and b, box plot of selected proteins based on the protein intensity before imputation
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
          ggpubr::stat_compare_means(comparisons = comparison_list, label="p.signif", paired = F,method="t.test") # This was replaced with a one-way ANOVA analysis, followed by a TukeyHSD test for post-hoc comparisons.
      
      list_merge.list[[protein]] = p1
  }
  
  merge.plot = ggpubr::ggarrange(plotlist = list_merge.list,ncol = 4,nrow = 1)
  pdf(paste0("write/Figure_S11.pdf"),width = 12,height = 2.5,onefile = T)
  print(merge.plot)
  dev.off()    
}  
## Obtain Pvalue from one-way ANOVA for selected 5 proteins
{
    rm(list=ls())
    Obj.list = readRDS(file="sp_dataset_03.rds")
    regions = c("Acinar","PanIN","PDAC")
    select_pid = c("Q7TPZ8","O09037","O88343","P07356","Q91YL7")
    meta.df = Obj.list$meta %>% dplyr::filter(Group %in% regions)
    meta.df$Regions = factor(meta.df$Group,levels = regions)
    meta.df = meta.df %>% arrange(Regions)
    input.df = Obj.list$filter[select_pid,meta.df$SampleID]
    
    list_merge.list= list()
    for(i in 1:dim(input.df)[1]){
        pid.df = input.df[i,] %>% t() %>% as.data.frame()
        names(pid.df) = "pid"
        pid.df$SampleId = row.names(pid.df)
        pid.df$Region = meta.df$Regions
        aov_model = aov(data=pid.df,pid~Region)
        p.value = summary(aov_model)[[1]]$`Pr(>F)`[1]

        tukey <- TukeyHSD(aov_model)
        print(select_pid[i])
        print(tukey) # Using this pvalue to replace
        print("===============")
    }
}
