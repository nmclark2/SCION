#####################################################################################
kmeans_clustering<-function(clustering_data, kmid) {
  
  set.seed(2020)
  #center and scale expression matrices
  normmatrix <- t(scale(t(clustering_data), scale=TRUE, center=TRUE)) 
  
  #test different clustering configurations, and use the silhouette index to pick the best one
  #cat("Calculating distance matrix\n")
  dist.mat <- dist(normmatrix)
  threshold <- NULL
  s.index <- NULL
  row <- 1
  kint <- floor(kmid*0.1) #how much to iterate k by
  #if kint < 8, change to 8 for a better starting number
  #min kmid is 20 so this works
  kint <- ifelse(kint<8,8,kint)
  for (k in seq(kmid-2*kint,kmid+2*kint,by=kint)){
    cat("Testing k=",k,"\n",sep="")
    #k-means clustering
    km.all <- suppressWarnings(kmeans(normmatrix, k,  nstart=5, iter.max=20))
    s <- silhouette(km.all$cluster,dist.mat)
    s.sum <- summary(s)
    threshold[row] <- k
    s.index[row] <- s.sum$avg.width
    row = row+1
  }
  
  #save the best configuration, which is the largest silhouette index 
  my.k <- threshold[s.index==max(s.index)]
  cat(paste("Choose k=",my.k,"\n",sep=""))
  km.all <- kmeans(normmatrix, my.k, nstart=50, iter.max=20)
  clusters <- as.data.frame(km.all$cluster)
  colnames(clusters) <- "clusters"
  results <- cbind(normmatrix,clusters)
  results <- as.data.frame(results)
  
  write.csv(results,'clusters.csv')
  return(results)
}
