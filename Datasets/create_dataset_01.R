#'@time 202403
#'@author Yuan
library(dplyr)

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
