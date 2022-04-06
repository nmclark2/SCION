#####################################################################################

dtw_clustering <- function(clustering_data,threshold)
{
#perform clustering using dtwclust
#check that the threshold is below 1
  # if it is not, set it to be 1
  if (threshold > 1){
    threshold = 1
  }
mydata = clustering_data
normmatrix <- (mydata[,]-rowMeans(mydata)*matrix(1,nrow=dim(mydata)[1],ncol=dim(mydata)[2]))/
  apply(mydata[,],1,sd)
normmatrix[is.na(rowSums(normmatrix)),] = rep(1,dim(normmatrix)[2])
#supervised clustering using dtw
clusters <- matrix(0,nrow=dim(normmatrix)[1],ncol=1)
ind=1
#for each gene, compare its profile to the rest of the genes
#any gene with a sufficiently matching profile is clustered together
for (j in 1:(dim(normmatrix)[1])){
  print(paste('Clustering gene ',toString(j),sep=""))
  #first check if gene j has already been clustered
  #if it has, we skip it
  #if it hasn't, we initialize the cluster #
  if (clusters[j]!=0){
    next
  }else{
    clusters[j]=ind
    ind = ind+1
  }
  for(k in (j+1):dim(normmatrix)[1]-1){
    #first check if gene k has already been clustered
    #if it has, we skip it
    if (clusters[k] != 0){
      next
    }
    #otherwise perform dtw on the two genes
    #DTW expects the time series as a vector
    ts1 = as.vector(t(normmatrix[j,]))
    ts2 = as.vector(t(normmatrix[k,]))
    
    #perform dtw
    dtwresults = dtw(ts1,ts2)
    
    #DTW gives us the indices that matched
    #if the time series matches perfectly, then x1=y1, x2=y2, etc.
    #we can test this by calculating the proportion of matches
    xinds = dtwresults$index1
    yinds = dtwresults$index2
    propmatches= sum(xinds==yinds)/length(xinds)
    
    #if propmatches>= threshold, cluster gene k with gene j
    #otherwise, move on
    if (propmatches>=threshold){
      clusters[k]=clusters[j]
    }
  }
}
#put together gene names, normalized expression, and clusters
results <- data.frame(normmatrix,clusters,row.names=row.names(mydata))
# #prepare data for plotting
# plotdata <- as.data.frame(t(results[,1:dim(normmatrix)[2]]))
# stacked <- stack(plotdata)
# stacked[,3] =  rep(colnames(results)[1:dim(normmatrix)[2]],dim(plotdata)[2])
# stacked[,4] = rep(clusters,each=dim(normmatrix)[2])
# colnames(stacked) <- c('Norm.Intensity','gene','group','cluster')
# #add means of each cluster as reference lines
# reflines = by(results[,1:dim(normmatrix)[2]], results$clusters, colMeans)
# reflinedata <- do.call(cbind, reflines)
# reflinedata <- as.data.frame(reflinedata)
# reflinestacked <- stack(reflinedata)
# reflinestacked[,2] <- rep(colnames(results)[1:dim(normmatrix)[2]],max(clusters))
# reflinestacked[,3] <- rep(1:max(clusters),each=dim(normmatrix)[2])
# colnames(reflinestacked) <- c('Norm.Intensity','group','cluster')
# #plot one cluster per file
# for (i in 1:max(results$clusters)){
#   g <- ggplot(data = stacked[stacked$cluster==i,],
#               mapping = aes(x=group, y=Norm.Intensity, colour = as.factor(cluster), group=gene)) +
#     geom_line() + theme(legend.position="none")+
#     scale_x_discrete(limits = unique(stacked$group), labels=colnames(mydata))
#   g <- g + geom_line(data = reflinestacked[reflinestacked$cluster==i,],
#                      aes(x = group, y = Norm.Intensity, group = cluster), colour='BLACK')
#   ggplotly(g)
#   ggsave(paste('Cluster Plots/cluster', i, '.png', sep=""),device="png",width=5, height=3,dpi=600)
# }
#finally, save the final data table
write.csv(results,'clusters.csv')
return(results)
}