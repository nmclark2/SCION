#Spatiotemporal Clustering for Integrative Omics Networks (SCION)

#target_genes_file: file containing list of target genes to use in network
#reg_genes_file: file containing list of regulator genes to use in network
#target_data_file: file containing table of data for targets (non-TFs) to use for inference
#reg_data_file: file containing table of data for regulators (TFs) to use for inference
#clustering_data_file: file containing table of data to use for both targets and regulators for clustering
#threshold: threshold to use for clustering (0<threshold<1)
#working_dir: link to the working directory you want to use

#For data tables, columns are samples are rows are genes. 
#Data are transposed before running the GENIE3 step. 

SCION <- function(target_genes_file,reg_genes_file,target_data_file,reg_data_file,clustering_data_file,threshold,working_dir){
  
  #check that all packages are installed, and load them
  for (package in c('dtwclust', 'plotly','randomForest')) {
    if (!require(package, character.only=T, quietly=T)) {
      install.packages(package)
      library(package, character.only=T)
    }else{
      library(package,character.only=T)
    }
  }
  
  #load functions and set working directory
  #make sure the functions are in your working directory!
  setwd(working_dir)
  
  #make folders we need
  dir.create('Cluster Plots')
  dir.create('Cluster networks')

  #read in tables
  target_genes = read.csv(target_genes_file,stringsAsFactors=FALSE)
  reg_genes = read.csv(reg_genes_file,stringsAsFactors=FALSE)
  target_data = read.csv(target_data_file,row.names=1)
  reg_data = read.csv(reg_data_file,row.names=1)
  clustering_data = read.csv(clustering_data_file,row.names=1)
  
  #get data for targets and regulators
  mytargetdata = target_data[row.names(target_data)%in%target_genes[,1],]
  myregdata = reg_data[row.names(reg_data)%in%reg_genes[,1],]
  
  #get clustering data
  myclustdata = clustering_data[row.names(clustering_data)%in%target_genes[,1],]
  
  #run clustering first
  clusterresults <- dtw_clustering(myclustdata,threshold)
  
  #infer a network using GENIE3 on each cluster
  for (i in 1:max(clusterresults$clusters)){
    mygenes = row.names(clusterresults[clusterresults$clusters==i,])
    clustertargetdata = mytargetdata[row.names(mytargetdata)%in%mygenes,]
    clusterregdata = myregdata[row.names(myregdata)%in%mygenes,]
    #if there's less than 2 TFs we skip the cluster
    #this is because there is no autoregulation
    if (dim(clusterregdata)[1]<=2){
      next
    }
    
    network = RS.Get.Weight.Matrix(t(clustertargetdata),t(clusterregdata))
    
    #trim matrix based on the ratio of TFs to targets - more TFs per target = more edges kept
    TFratio = (dim(clusterregdata)[1])/(dim(clustertargetdata)[1])
    if (TFratio < 0.25){
      edgestokeep = floor(0.5*dim(clustertargetdata)[1])
    } else if (TFratio >= 0.25 & TFratio < 0.5){
      edgestokeep = floor(1.5*dim(clustertargetdata)[1])
    } else{
      edgestokeep = floor(2*dim(clustertargetdata)[1])
    }
    
    #get the weight threshold that corresponds to this number of edges to keep
    allweights = stack(data.frame(network))
    allweights = allweights[order(allweights[,1],decreasing=TRUE),1]
    #its possible the number of edges to keep has a threshold of 0 or NA
    #if this is the case we all edges with weights >0
    if(is.na(allweights[edgestokeep]) | allweights[edgestokeep]==0){
      weightthreshold = min(na.omit(allweights[allweights>0]))
    }else{
      weightthreshold = allweights[edgestokeep]
    }
    
    #now make a new network where we eliminate all the low confidence edges
    trimmednet = data.frame(ifelse(network<weightthreshold,NaN,network))
    
    #translate the trimmed network into a table we can import into cytoscape
    networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(),
                              stringsAsFactors=FALSE)
    row=1
    for (j in 1:dim(trimmednet)[1]){
      for (k in 1:dim(trimmednet)[2]){
        #skip NaNs as these have no edge
        if (is.na(trimmednet[j,k])){
          next
        }else{
          networktable[row,] = cbind(colnames(trimmednet)[k],"regulates",rownames(trimmednet)[j],trimmednet[j,k])
          row = row+1
        }
      }
    }
    
    #if there are . in the networktable, replace them with - (for phospho)
    if (length(grepl('\\.',networktable$Regulator))>0){
      networktable$Regulator <- sapply(networktable$Regulator, function(x) gsub("\\.", "-", x))
    }
    
    #write table to file
    write.table(networktable,paste('Cluster networks/network_cluster_',i,'.txt',sep=""),row.names=FALSE,quote=FALSE,sep='\t')
    
    #get the hub gene (highest outdegree) and save it to connect the clusters later
    #if there is a tie, we save both hubs
    myregs = unique(networktable$Regulator)
    for (j in 1:length(myregs)){
      numedges = sum(networktable$Regulator%in%myregs[j])
      if (j==1){
        hubedges = numedges
        hub = myregs[j]
      }else if (numedges>hubedges){
        hubedges = numedges
        hub = myregs[j]
      }else if (numedges==hubedges){
        hub = rbind(hub,myregs[j])
      }
    }
    if (!exists("myhubs")){
      myhubs = hub
    } else{
      myhubs = rbind(myhubs,hub)
    }
  }
  #infer network connecting the hubs
  #as before we can only do this if there is more than 1 hub
  if (exists("myhubs") && length(myhubs)>2){
    hubtargetdata = mytargetdata[row.names(mytargetdata)%in%myhubs,]
    hubregdata = myregdata[row.names(myregdata)%in%myhubs,]
    network = RS.Get.Weight.Matrix(t(hubtargetdata),t(hubregdata))
    
    #trim matrix based on the ratio of TFs to targets - more TFs per target = more edges kept
    edgestokeep = floor(4*dim(hubtargetdata)[1])
    
    #get the weight threshold that corresponds to this number of edges to keep
    allweights = stack(data.frame(network))
    allweights = allweights[order(allweights[,1],decreasing=TRUE),1]
    weightthreshold = allweights[edgestokeep]
    
    #now make a new network where we eliminate all the low confidence edges
    trimmednet = data.frame(ifelse(network<weightthreshold,NaN,network))
    
    #translate the trimmed network into a table we can import into cytoscape
    networktable = data.frame(Regulator=character(), Interaction=character(), Target=character(), Weight=double(),
                              stringsAsFactors=FALSE)
    row=1
    for (j in 1:dim(trimmednet)[1]){
      for (k in 1:dim(trimmednet)[2]){
        #skip NaNs as these have no edge
        if (is.na(trimmednet[j,k])){
          next
        }else{
          networktable[row,] = cbind(colnames(trimmednet)[k],"regulates",rownames(trimmednet)[j],trimmednet[j,k])
          row = row+1
        }
      }
    }
    
    #if there are . in the networktable, replace them with - (for phospho)
    if (length(grepl('\\.',networktable$Regulator))>0){
      networktable$Regulator <- sapply(networktable$Regulator, function(x) gsub("\\.", "-", x))
    }
    
    #write table to file
    write.table(networktable,'Cluster networks/network_hub.txt',row.names=FALSE,quote=FALSE,sep='\t')
  }
  
  #read in all the tables, combine, and save as the final file to import into cytoscape
  #get all of the folders for our samples
  myfiles = list.files(path='Cluster networks')
  for (i in 1:length(myfiles)){
    if (i==1){
      finalnetwork = read.table(paste('Cluster networks/',myfiles[i],sep=""),header=TRUE)
    }else{
      finalnetwork = rbind(finalnetwork, read.table(paste('Cluster networks/',myfiles[i],sep=""),header=TRUE))
    }
  }
  write.table(finalnetwork,'FINAL_NETWORK.txt',row.names=FALSE,quote=FALSE,sep='\t')
  message('Finished!')
}
