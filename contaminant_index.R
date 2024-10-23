Input_dir<-list.files(paste0(getwd(), "/wkdir"))
Input_dir

#import_files
library(dplyr)
library(readr)
library(stringr)

#######import csv files#######
#
for (j in 1:length(Input_dir)) {
  
  ###
  My_path<-list.files(paste0(getwd(), "/wkdir/",Input_dir[j]),recursive = TRUE, pattern = "(.*)csv$")
  My_path<-paste0(Input_dir[j],"/",My_path)
  My_path
  ##the path of csv files
  csv_path<-paste0(getwd(), "/wkdir/",My_path)
  #store csv file into a list
  csv_list<-lapply(csv_path, read_csv)
  names(csv_list)
  #
  csv_names<-str_remove_all(list.files(paste0(getwd(), "/wkdir/",Input_dir[j]),recursive = TRUE, pattern = "(.*)csv$"), 
                            ".csv")
  #
  names(csv_list)<-csv_names
  
  ###build a dataframe with filenames, RT and contamination index##
  #length(csv_list)
  contaminant_index_inten<-data.frame(file_name=names(csv_list)
                                      ,rt_min=array(dim = length(csv_list)),
                                      tot_inten=array(dim = length(csv_list)),
                                      con_index=array(dim = length(csv_list))
  )
  for (i in 1:length(csv_list)) {
    Temp<-csv_list[[i]]
    Temp<-Temp%>%
      select(rt_values_min, mz_values, mobility_values, corrected_intensity_values)
    Temp$count<-ifelse(0.00076923*Temp$mz_values-Temp$mobility_values+0.56923>0, 1, c(-1))
    Temp$inten<-Temp$corrected_intensity_values*Temp$count
    
    ##contaminant index
    #n_con<-Temp%>%
    #filter(count==c(-1))
    #n_con<-abs(sum(n_con$count))
    #n_sig<-Temp%>%
    #filter(count==c(1))
    #n_sig<-abs(sum(n_sig$count))
    #contaminant_index<-(n_con/n_sig)*100
    ##
    n_con<-Temp%>%
      filter(count==c(-1))
    n_con<-abs(sum(n_con$inten))
    n_sig<-Temp%>%
      filter(count==c(1))
    n_sig<-abs(sum(n_sig$inten))
    contaminant_index_inten[i,"rt_min"]<-as.numeric(unique(csv_list[[i]]["rt_values_min"])) 
    contaminant_index_inten[i,"con_index"]<-(n_con/n_sig)*100
    contaminant_index_inten[i,"tot_inten"]<-as.numeric(log2(sum(csv_list[[i]]["rt_values_min"]))) 
  }
  
  #plot
  plot(contaminant_index_inten$rt_min, contaminant_index_inten$con_index,
       xlim =c(0,80),
       ylim=c(0,200)
  )
  
  plot(contaminant_index_inten$rt_min, contaminant_index_inten$tot_inten,
       xlim =c(0,80),
       ylim=c(5,30)
  )
  
  #####organize and export########
  #
  ###export results in output folder and set the file name as the [experiment name].csv
  if (!dir.exists(paste0(getwd(),"/output")))
    dir.create(paste0(getwd(),"/output")) #create output folder
  output<-paste0(paste0(getwd(),"/output/"),
                 Input_dir[j],".csv") #set the file name as the [experiment name].csv
  Temp<-contaminant_index_inten
  write.csv(Temp, output)
}





