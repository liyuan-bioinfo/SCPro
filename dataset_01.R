#'@time 202403
#'@author Yuan

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(enrichplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)
setwd("") 

# -------------------------------------------------------------------------
#                            I-Basic pre-treat                            #
# -------------------------------------------------------------------------
{
    # Desc: Create .RDS and save files from Perseus output, including Quantified, Filter and Impute
    # Output: .RDS
    {
        rm(list=ls())
        # prepare meta files
        meta_df = read.delim("meta.txt",header=T,check.names=FALSE) 
        ct3_order = c("Acinar","Lymph","Tumor")
        meta_df$CellType = factor(meta_df$CellType, levels=ct3_order)
        meta_df = meta_df %>% arrange(CellType)

        # prepare Quantified file
        quantified_df = read.delim("Quantified_2917.txt",header=T,check.names=FALSE) #6406 * 69
        quantified_df = quantified_df[,meta_df$SampleID]#re-order
        names(quantified_df) = meta_df$SampleID #re-name

        # prepare Impute file
        impute_df = read.delim("Impute_2608.txt",header=T,check.names=FALSE) #5900 * 69
        row.names(impute_df) = impute_df$Protein
        impute_df = impute_df[,meta_df$SampleID]
        names(impute_df) = meta_df$SampleID

        # prepare Impute_mean file
        impute_mean_df = read.delim("Impute_mean_2608.txt",header=T,check.names=FALSE) #5900 * 14
        row.names(impute_mean_df) = impute_mean_df$Protein
        impute_mean_df = impute_mean_df[,ct3_order]

        # save  
        Obj_list = list() # save RDS
                
        Obj_list$quantified_df = quantified_df
        Obj_list$impute_df = impute_df
        Obj_list$impute_mean_df = impute_mean_df
        Obj_list$meta_df = meta_df
        Obj_list$ct3_order = ct3_order
        
        saveRDS(Obj_list, file="dataset_01.rds")       
    }

    # Desc: ANOVA analysis of celltype groups, including 14 cell types and 4 cell lineages
    # Input: Obj_list$impute_df; Obj_list$impute_mean_df
    # Output: Obj_list$aov_ct14_df
    # for Fig.5e,f;  Fig.6c,d, Extended Data Fig. 10 and Fig. 12
    {
        rm(list=ls())
        Obj_list = readRDS(file="dataset_01.rds") 

        meta_df = Obj_list$meta_df        
        data_df = Obj_list$impute_df
        data_mean_df = Obj_list$impute_mean_df
        
        ct_num = 3        
        protein_num = dim(data_df)[1]
        ## log2 FC, to find enriched cell-types with mean abundance
        dep_df = data.frame()
        for(i in 1:ct_num){
            temp_df = data_mean_df[,i] - data_mean_df[,1:ct_num]
            temp_pid = names(which(rowSums(temp_df>log2(1)) == (ct_num-1)))
            temp_df2 = data.frame(pids=temp_pid)
            temp_df2$CellType = names(data_mean_df)[i]
            temp_df2$log2FC = apply(temp_df[temp_pid,-i],1,mean)        
            dep_df = rbind(dep_df,temp_df2)            
        }

        ## pvalue, with each sample
        row.names(dep_df) = dep_df$pids
        dep_df = dep_df[row.names(data_df),] # re-order

        Pvalue = c()
        for(i in 1:protein_num){
            pid_df = data_df[i,] %>% t() %>% as.data.frame()
            names(pid_df) = "pid"
            pid_df$SampleId = row.names(pid_df)
            pid_df$CellType = meta_df$CellType
            model = aov(data=pid_df,pid~CellType)
            # temp_pvalue = summary(model)[[1]]$`Pr(>F)`[1]
            sig_df = TukeyHSD(model)$CellType #which ct is enriched
            enriched_ct = dep_df[i,"CellType"]
            enriched_ct_sig_df = sig_df[grep(row.names(sig_df),pattern=paste0("-",enriched_ct,"$|","^",enriched_ct,"-")),] #whethor this ct is sig. with other cts
            temp_pvalue = RecordTest::fisher.method(as.numeric(enriched_ct_sig_df[,"p adj"]))$"p.value"[[1]]  
            # temp_pvalue = max(enriched_ct_sig_df[,"p adj"])                                    
            Pvalue = c(Pvalue,temp_pvalue)            
        }        

        dep_df$pvalue = Pvalue
        dep_df$fdr = p.adjust(Pvalue,method = "BH")
        # write.csv(dep_df,"write/Stat/dep_ANOVA.csv")
        # head(dep_df)

        temp_dep_df = dep_df %>% filter(log2FC>=1 & fdr<0.05)
        table(temp_dep_df$CellType)
        Obj_list$aov_ct3_df = dep_df

        saveRDS(Obj_list, file="dataset_01.rds")
    }        

}
# -------------------------------------------------------------------------
#                            II-Quality Control                           #
# -------------------------------------------------------------------------
{
    # PCA analysis
    # Input: Obj_list$quantified_df; Obj_list$meta_df;
    # Fig.2i
    {
        rm(list=ls())
        Obj_list = readRDS("dataset_01.rds")
        
        data_df = Obj_list$quantified_df
        meta_df = Obj_list$meta_df
        
        ## PCA analysis
        data_df = log2(data_df+1) %>% t() %>% as.data.frame() 

        pca_result <- prcomp(data_df,scale=T,center = T)
        plot_df = as.data.frame(pca_result$x)
        summ1 <- summary(pca_result)
        xlab1 <- paste0("PC1 (",round(summ1$importance[2,1]*100,2),"%)")
        ylab1 <- paste0("PC2 (",round(summ1$importance[2,2]*100,2),"%)")

        p1=ggplot(data = plot_df,aes(x = PC1,y = PC2,color = meta_df$CellType))+
            stat_ellipse(aes(fill = meta_df$CellType),
                            type = "norm",geom = "polygon",alpha = 0.25, color = NA)+ 
            geom_point(size = 2)+
            labs(x = xlab1,y = ylab1)+ guides(fill = "none")+ theme_bw()+            
            theme(plot.background = element_blank(),legend.background = element_blank(),
                panel.background = element_blank(),panel.grid = element_blank(),
                axis.text = element_text(size = 12),axis.title = element_text(size = 16),
                legend.text = element_text(size = 16))

        ## save plot
        pdf(file=paste0("write/Quantified/Fig_2i_PCA_",Sys.Date(),".pdf"),
            width = 7,height=5)
        print(p1)
        dev.off()    

        ## save table        
        write.csv(plot_df, file = paste0("write/Quantified/Fig_2i_PCA_PCA_",Sys.Date(),".csv"))
    }


    # Correlation analysis
    # Input: Obj_list$quantified_df; Obj_list$meta_df;
    # ext Fig.2b
    {
        rm(list=ls())
        Obj_list = readRDS("dataset_01.rds")
        
        data_df = Obj_list$impute_df
        meta_df = Obj_list$meta_df
                
        # data_df = log2(data_df+1)
        plot_df = cor(data_df[])
        
        ## save table        
        write.csv(plot_df, file = paste0("write/Quantified/sFig_2b_Corr_",Sys.Date(),".csv"))

        # search the max value
        row_max = c()
        for(i in 1:dim(plot_df)[1]){
            temp_row = plot_df[i,]
            row_max = c(row_max,max(temp_row[which(temp_row!=1)]))            
            
        }

        # replace "1" with max corr
        row_max = max(row_max)
        for(i in 1:dim(plot_df)[1]){
            temp_row = plot_df[i,]            
            plot_df[i,which(temp_row==1)]=row_max            
        }    

        p1 = pheatmap(plot_df,cluster_rows = T,cluster_cols = T,width = 7,height = 6,#cellwidth = 3,cellheight = 3,
                        show_rownames = T,show_colnames = T,border_color = NA,
                        annotation_names_row = F,annotation_names_col = F,                        
                        annotation_legend = F,fontsize_row = 10,fontsize_col = 10,
                        color = colorRampPalette(c("blue","white","red"))(100),silent = TRUE)

        # add legends
        llim = round(min(plot_df),2)
        ulim = round(max(plot_df),2)

        # Fig lengends           
        df <- data.frame(matrix(nrow = 10, ncol = 2))
        df[] <- rnorm(20)
        p2=ggplot(data = df, aes(x = X1, y = X2, fill = X2)) + geom_tile() + theme_void()+
            scale_fill_gradientn(colors=colorRampPalette(c("blue","white","red"))(100),
            limits = c(llim, ulim), breaks = round(seq(llim, ulim, length.out=3),2),guide = guide_colorbar(barwidth = 5, barheight = 1, 
                                                                                            ticks = TRUE, ticks.colour="black",ticks.linewidth=1/.pt,direction="horizontal",#"vertical",                                                                                            
                                                                                            label.theme = element_text(size = 8),
                                                                                            draw.ulim = TRUE, 
                                                                                            draw.llim = TRUE, 
                                                                                            draw.separator = TRUE,
                                                                                            
                                                                                            title.position = "top"
                                                                                            )
            )     

        ## save plot
        pdf(file=paste0("write/Quantified/sFig_2b_Corr_",Sys.Date(),".pdf"),
            width = 7,height=5)
        print(p1)
        print(p2)
        dev.off()    


    }    
}

# -------------------------------------------------------------------------
#                            III-Stat Analysis                            #
# -------------------------------------------------------------------------
{   
    # Desc: Pheatmap of sig.protein panel of 3 cts
    # Input: Obj_list$aov_ct3_df;
    # Output: Pheatmap[Fig.2j]; Related table
    {
        rm(list=ls())
        Obj_list = readRDS(file="dataset_01.rds") 

        meta_df = Obj_list$meta_df                
        data_df = Obj_list$impute_df[,meta_df$SampleID]

        dep_df = Obj_list$aov_ct3_df %>% dplyr::filter(fdr<0.05 & log2FC>1) # 596
        
        dep_df$CellType = factor(dep_df$CellType, levels = Obj_list$ct3_order)
        dep_df = dep_df %>% dplyr::arrange(CellType)
        
        plot_df = data_df[dep_df$pids,]                                        
                
        # 使用apply函数对data.frame的每一行应用缩放函数  
        df_scaled_rows <- t(apply(plot_df, 1, scale))  
        
        # 计算99%分位数
        ulim <- quantile(df_scaled_rows, 0.99)
        llim <- quantile(df_scaled_rows, 0.01)

        # 将超过99%分位数的值替换为99%分位数
        df_scaled_rows <- ifelse(df_scaled_rows > ulim, ulim, df_scaled_rows)
        df_scaled_rows <- ifelse(df_scaled_rows < llim, llim, df_scaled_rows)

        # 将结果转换回data.frame格式  
        df_scaled_rows <- as.data.frame(df_scaled_rows)
        names(df_scaled_rows) = names(plot_df)

        length(which(df_scaled_rows<0))
        length(which(df_scaled_rows>0))        

        # 创建颜色映射
        my_palette <- c(
            c(colorRampPalette(c("blue", "white"))(51)),
            c(colorRampPalette(c("white", "red"))(40))
        )
        # 将0设置为断点
        my_breaks <- unique(c(seq(min(df_scaled_rows), 0, length.out = 51), seq(0, max(df_scaled_rows), length.out = 40)))

        p1=pheatmap::pheatmap(df_scaled_rows,cluster_rows = F,show_rownames = F,cluster_cols = F,
                scale = "none",border_color = NA,legend = T,cellwidth = 5,cellheight=0.3,              
                color = my_palette,breaks=my_breaks,                
                show_colnames = F,legend_labels = F,annotation_legend = F,labels_col = NA,labels_row = NA,
                annotation_names_col = FALSE,silent = TRUE
        )

        # Fig lengends           
        my_palette <- c(
            c(colorRampPalette(c("blue", "white"))(50)),
            c(colorRampPalette(c("white", "red"))(50))
        )        
        df <- data.frame(matrix(nrow = 10, ncol = 2))
        df[] <- rnorm(20)
        p2=ggplot(data = df, aes(x = X1, y = X2, fill = X2)) + geom_tile() + theme_void()+#
            scale_fill_gradientn(colors=my_palette, limits = c(llim, ulim), breaks = c(-1.5,0,1.5),guide = guide_colorbar(barwidth = 5, barheight = 1, 
                                                                                            ticks = TRUE, ticks.colour="black",ticks.linewidth=1/.pt,direction="horizontal",#"vertical",                                                                                            
                                                                                            label.theme = element_text(size = 8),
                                                                                            draw.ulim = TRUE, 
                                                                                            draw.llim = TRUE, 
                                                                                            draw.separator = TRUE,
                                                                                            
                                                                                            title.position = "top"
                                                                                            )        
            )

        pdf(file=paste0("write/Stat/Fig2_j_Pheatmap_",Sys.Date(),".pdf"))
        print(p1)
        print(p2)
        dev.off()                                                                                            

    }
}

# -------------------------------------------------------------------------
#                          IV-Functional Enrichment analysis              #
# -------------------------------------------------------------------------
{
    # Desc: GO enrichment of cell-type enriched proteins for four major cell-lineages
    # Input: 
    # Output: 
    {
        rm(list=ls())        
        Obj_list = readRDS(file="dataset_01.rds") 
        meta_df = Obj_list$meta_df                

        dep_df = Obj_list$aov_ct3_df %>% dplyr::filter(fdr<0.05 & log2FC>1) # 596        
        dep_df$pid = gsub(dep_df$pids,pattern=";.*",replacement="")        

        tran_df = bitr(dep_df$pids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
        plot_df = dep_df %>% merge(tran_df, by.y="UNIPROT", by.x="pid",all.x=F,all.y=F)    
        plot_df$CellType = factor(plot_df$CellType, levels=Obj_list$ct3_order)    
        plot_df = plot_df %>% dplyr::arrange(CellType)

        formula_res <- compareCluster(data = plot_df,ENTREZID~CellType,fun="enrichGO", 
                                    OrgDb = org.Mm.eg.db,
                                    ont = "BP", pAdjustMethod = "BH",
                                    pvalueCutoff = 1,qvalueCutoff = 1,readable=TRUE
                                    
        )
        formula_res_cutoff = formula_res
        formula_res_cutoff@compareClusterResult = formula_res@compareClusterResult[formula_res@compareClusterResult$p.adjust<=0.05,]

        write.csv(x = formula_res@compareClusterResult, file=paste0("write/Functional/Fig2_j_GOBP_",Sys.Date(),".csv"))  
                   
    }

}
