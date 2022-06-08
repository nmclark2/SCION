#NMS Calculation Code
#March 2022
#This assumes the motifs are in a text file output from the NetMatchStar App in Cytoscape. The code may need to be modified for other outputs.
#The node table from the network is also needed to run this code.

#read in the node table
nodes <- read.csv('node-table.csv',stringsAsFactors=FALSE)
nodes <- nodes$shared.name

#for each of the motifs, count the number of motifs for each gene and save in a column
nms <- data.frame(nodes)
motifs <- list.files(pattern="*.txt")
for(i in 1:length(motifs)){
  print(i)
  my.motif <- read.table(motifs[i],header=FALSE,skip=4,sep='\t')
  my.motif.genes = unlist(strsplit(my.motif[,5],', '))
  motif.num=NULL
  for(j in 1:dim(nms)[1]){
    motif.num <- rbind(motif.num,cbind(nodes[j],sum(my.motif.genes%in%nodes[j])))
  }
  nms <- cbind(nms,as.numeric(motif.num[,2]))
  colnames(nms)[i+1] = make.names(unlist(strsplit(motifs[i],".txt")))
}

#calculate the NMS
#first, normalize the values for each column to a 0 to 1 scale
norm.cols <- NULL
for(i in 1:length(motifs)){
  my.col <- nms[,i+1]
  my.col.norm <- (my.col-min(my.col))/(max(my.col)-min(my.col))
  norm.cols <- cbind(norm.cols,my.col.norm)
}
colnames(norm.cols) <- paste("norm.",colnames(nms)[2:(length(motifs)+1)],sep="")

#next, calculate the multiplier, which is the number of motifs each gene is in
multiplier <- rowSums(norm.cols>0)

#finally, calculate the NMS, which is the sum of the normalized motifs times the multiplier
nms.score <- rowSums(norm.cols)*multiplier

#combine everything and save
results <- cbind(nms,norm.cols,multiplier,nms.score)
write.csv(results,'NMS.csv',row.names=FALSE)

