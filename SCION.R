#Version 2.0 beta
#July 2020

SCION <- function(target_genes_file,reg_genes_file,target_data_file,reg_data_file,
                  is_clustering,clustering_data_file,threshold,clusters_file,hubs,working_dir){
  
  #set working directory which contains all your input files
  setwd(working_dir)

  #if clustering, create folders
  if (!grepl("None",is_clustering)){
    #make folders we need
    dir.create('Cluster Plots')
    dir.create('Cluster networks')
  }

  #read in tables
  target_genes = read.csv(target_genes_file,stringsAsFactors=FALSE)
  reg_genes = read.csv(reg_genes_file,stringsAsFactors=FALSE)
  target_data = read.csv(target_data_file,row.names=1)
  reg_data = read.csv(reg_data_file,row.names=1)
  
  #get data for targets and regulators
  mytargetdata = target_data[row.names(target_data)%in%target_genes[,1],]
  myregdata = reg_data[row.names(reg_data)%in%reg_genes[,1],]
  
  #read clustering data if needed
  if (grepl("Temporal",is_clustering)){
    clustering_data = read.csv(clustering_data_file,row.names=1)
    myclustdata = clustering_data[row.names(clustering_data)%in%target_genes[,1],]
  }
  
  #run clustering first
  #if inputting a clusters file, read it in here
  if (is_clustering=="Temporal"){
    clusterresults <- dtw_clustering(myclustdata,threshold)
  }else if (is_clustering=="Non-Temporal"){
    clusterresults <- ica_clustering(myclustdata,threshold)
  }else if (is_clustering=="Upload"){
    clusterresults <- read.csv(clusters_file,row.names=1)
  }
  
  #infer a network using GENIE3 on each cluster
  #if not clustering, just infer one network
  if (is_clustering=="None"){
    mygenes = target_genes
    clustertargetdata = mytargetdata
    clusterregdata = myregdata
    
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
      print(j)
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
    write.table(networktable,'FINAL_NETWORK.txt',row.names=FALSE,quote=FALSE,sep='\t')
    message('Finished!')
  }
  #otherwise, infer 1 network for each cluster
  else{
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
    #if desired, infer network connecting the hubs
    #as before we can only do this if there is more than 1 hub
    if (hubs=="Yes"){
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
}
