#example script to pre-process data tables for SC-ION
library(pacman)
p_load(openxlsx)
p_load(stringr)
p_load(ExPosition)

#file paths
#target matrix - data we want to use for our targets (typically transcript)
target.str <- 'JA_FPKM.csv'

#cluster matrix - data we want to use for clustering
cluster.str <- 'JA_FPKM_means.csv'

#list of regulators - typically transcription factors
TF.str <- 'Arabidopsis_TFs_AGRIS.xlsx'

#list of candidates/targets
candidates <- 'SC-ION_candidates_updated.txt'

#read in files
target.mat <- read.csv(target.str,row.names=1)
cluster.mat <- read.csv(cluster.str,row.names=1)
#need to sort cluster.mat by column names
cluster.mat <- cluster.mat[,order(colnames(cluster.mat))]
TF <- read.xlsx(TF.str)
candidates <- read.table(candidates)

#process target.mat first
#here, we only take the 0 and 2hr timepoints
my.target.mat <- target.mat[,grepl("hr0_",colnames(target.mat)) | grepl("hr2_",colnames(target.mat))]
#row normalize and remove missing values
my.target.mat <- rowNorms(my.target.mat, type="z")
my.target.mat <- na.omit(my.target.mat)
#save
write.csv(my.target.mat,'target_mat_RNA.csv',na="0")
#save target list
my.targets <- row.names(my.target.mat)[row.names(my.target.mat)%in%candidates[,1]]
write.csv(my.targets,'target_list_RNA.csv',row.names=FALSE)

#function to process regulator matrix
#accounts for whether input data are phospho or protein
process_reg_mat <- function(reg.str, phospho, targets, cluster.mat, TF, candidates){
  
  #read in regulator matrix
  reg.mat <- read.csv(reg.str,row.names=1)
  
  #next process the regulator matrix
  #different depending on protein vs phospho
  if(phospho){
    #merge multiplicities
    mat <- aggregate(reg.mat, list(reg.mat$original.id), function(x) mean(x, na.rm=T))
    my.data <- mat[,grepl("JA",colnames(mat)) | colnames(mat)%in%"original.id"]
    #get site information for each ID
    reg.mat$original.id <- as.numeric(reg.mat$original.id)
    reg.mat <- reg.mat[order(reg.mat$original.id),]
    my.ids <- reg.mat$original.id[!duplicated(reg.mat$original.id)]
    my.proteins <- reg.mat$Protein[!duplicated(reg.mat$original.id)]
    my.sites <- reg.mat$Positions.within.proteins[!duplicated(reg.mat$original.id)]
    my.names <- NULL
    for(i in 1:length(my.proteins)){
      #remove isoform information
      my.pro <- unlist(strsplit(my.proteins[i],'\\.'))[1]
      #get site info
      my.site <- unlist(strsplit(my.sites[i],';'))[1]
      #paste and save
      my.id <- paste(my.pro,'.p',my.site,sep="")
      my.results <- data.frame(my.id,my.ids[i])
      my.names <- rbind(my.names,my.results)
    }
    colnames(my.names) <- c('Protein.Site','original.id')
    #merge
    final.data <- merge(my.names,my.data,by="original.id",sort=FALSE)
    #collapse duplicate sites
    if(sum(duplicated(final.data$Protein.Site)) > 0){
      mat <- aggregate(final.data, list(final.data$Protein.Site), function(x) mean(x, na.rm=T))
      mypros <- mat$Group.1
      mat <- mat[, -c(1:3)]
      final.data <- mat
      row.names(final.data) <- mypros
    }else{
      row.names(final.data) <- final.data$Majority.protein.IDs
      final.data <- final.data[,-c(1:2)]
    }
    #row normalize and save
    final.data <- rowNorms(final.data, type="z")
    write.csv(final.data,'reg_mat_phospho.csv',na="0")
    
    #save the regulator list
    TF.list <- TF[TF[,1]%in%candidates[,1],1]
    my.genes <- unlist(strsplit(row.names(final.data),'\\.'))
    my.genes <- my.genes[grepl('AT',my.genes)]
    TF.final <- row.names(final.data)[my.genes%in%TF.list]
    write.csv(TF.final,'reg_list_phospho.csv',row.names=FALSE)
  }else{
    #use majority protein
    my.proteins <- data.frame(reg.mat$Majority.protein.IDs,as.numeric(reg.mat$id))
    colnames(my.proteins) <- c('Majority.protein.IDs','id')
    #if multiple IDs, stack
    my.proteins.stack <- NULL
    for(i in 1:dim(my.proteins)[1]){
      my.id <- my.proteins$Majority.protein.IDs[i]
      if(grepl(';',my.id)){
        my.ids <- unlist(strsplit(my.id,';'))
        my.ids <- my.ids[grepl('AT',my.ids)]
        #remove isoform information
        my.ids <- unlist(strsplit(my.ids,'\\.'))
        my.ids <- my.ids[grepl("AT",my.ids)]
        #remove duplicates
        my.ids <- unique(my.ids)
        my.results <- data.frame(my.ids,rep(my.proteins$id[i],length(my.ids)))
        colnames(my.results) <- c('Majority.protein.IDs','id')
        my.proteins.stack <- rbind(my.proteins.stack,my.results)
      }else{
        #remove isoform information
        my.id <- unlist(strsplit(my.id,'\\.'))
        my.id <- my.id[grepl("AT",my.id)]
        my.results <- data.frame(my.id,my.proteins$id[i])
        colnames(my.results) <- c('Majority.protein.IDs','id')
        my.proteins.stack <- rbind(my.proteins.stack,my.results)
      }
    }
    #get expression data
    #for this experiment, take JA only
    my.data <- reg.mat[,grepl("JA",colnames(reg.mat)) | colnames(reg.mat)%in%"id"]
    my.data$id <- as.numeric(my.data$id)
    final.data <- merge(my.proteins.stack,my.data,by="id",sort=FALSE)
    #aggregate expression for duplicated gene symbols
    if(sum(duplicated(final.data$Majority.protein.IDs)) > 0){
      mat <- aggregate(final.data, list(final.data$Majority.protein.IDs), function(x) mean(x, na.rm=T))
      mypros <- mat$Group.1
      mat <- mat[, -c(1:3)]
      final.data <- mat
      row.names(final.data) <- mypros
    }else{
      row.names(final.data) <- final.data$Majority.protein.IDs
      final.data <- final.data[,-c(1:2)]
    }
    #row normalize and save
    final.data <- rowNorms(final.data, type="z")
    write.csv(final.data,'reg_mat_protein.csv',na="0")
    #save the regulator list
    TF.list <- TF[TF[,1]%in%candidates[,1],1]
    TF.final <- TF.list[TF.list%in%row.names(final.data)]
    write.csv(TF.final,'reg_list_protein.csv',row.names=FALSE)
  }
  
  #make the clustering matrix
  #the clustering matrix is RNA
  #for phospho, we have to add on the phosphosites as additional rows with their corresponding clustering data
  #otherwise, they will not be in the clustering matrix
  if(phospho){
    #get data for targets
    cluster.data <- cluster.mat[row.names(cluster.mat)%in%my.targets,]
    #for each TF, get the data for that TF and append to the cluster matrix
    for (i in 1:length(TF.final)){
      my.TF <- unlist(strsplit(TF.final[i],'\\.'))[1]
      cluster.data <- rbind(cluster.data,cluster.data[row.names(cluster.data)%in%my.TF,])
      rownames(cluster.data)[dim(cluster.data)[1]]=TF.final[i]
      #row normalize and save
      cluster.data <- rowNorms(cluster.data, type="z")
      write.csv(cluster.data,'cluster_mat_phospho.csv',na="0")
    }
  }else{
    #for the protein, we just match 1:1, row normalize, and save
    cluster.data <- cluster.mat[row.names(cluster.mat)%in%my.targets | row.names(cluster.mat)%in%TF.final,]
    cluster.data <- rowNorms(cluster.data, type="z")
    write.csv(cluster.data,'cluster_mat_protein.csv',na="0")
  }
}

#process proteome data
process_reg_mat('Normalized_values_protein.csv',FALSE,my.targets, cluster.mat,TF, candidates)

#process phospho data
process_reg_mat('Normalized_values_phospho.csv',TRUE, my.targets, cluster.mat,TF, candidates)