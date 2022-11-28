#Original author Mitch Elmore
#Modified by Natalie Clark for SCION

#### ICAclust ####
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0181195

#smaller values of k = more and tighter clusters

#####################################################################################
ica_clustering<-function(clustering_data, k) {
  
  normmatrix <- t(scale(t(clustering_data), scale=TRUE, center=TRUE)) #center and scale expression matrix
  
  #for ICA, smaller k = tighter clusters, so we need to convert the threshold
  #k = 1.1-k
  
  #set.seed(2020)
  X.ICA<-fastICA(normmatrix,n.comp=ncol(normmatrix), alg.typ = "parallel", fun = "logcosh", alpha = 1.0,
                 method = "C", row.norm = FALSE, maxit = 5000, tol = 1e-03, verbose = TRUE)
  hc.ICA<-hclust(dist(X.ICA$S), method = "ward.D", members=NULL)
  mojena=mean(hc.ICA$height)+k*sd(hc.ICA$height)
  cluster_num=length(hc.ICA$height[hc.ICA$height>mojena]) + 1
  clusters<-cutree(hc.ICA, k=cluster_num)
  results <- cbind(normmatrix,clusters)
  
  #convert matrix to dataframe
  results <- as.data.frame(results)
  
  #prepare data for plotting
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
  #     scale_x_discrete(limits = unique(stacked$group), labels=colnames(normmatrix))
  #   g <- g + geom_line(data = reflinestacked[reflinestacked$cluster==i,],
  #                      aes(x = group, y = Norm.Intensity, group = cluster), colour='BLACK')
  #   ggplotly(g)
  #   ggsave(paste('Cluster Plots/cluster', i, '.png', sep=""),device="png",width=5, height=3,dpi=600)
  # }
  #finally, save the final data table
  write.csv(results,'clusters.csv')
  return(results)
}
