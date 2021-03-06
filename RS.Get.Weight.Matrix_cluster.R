# target.matrix and input.matrix (TFs) have samples as rows and genes as columns
RS.Get.Weight.Matrix_cluster<- function(target.matrix, input.matrix, K="sqrt", nb.trees=10000, importance.measure="%IncMSE", seed=NULL, trace=TRUE, ...)  
{
  require(parallel)
  # set random number generator seed if seed is given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # to be nice, report when parameter importance.measure is not correctly spelled
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }
  
  # normalize expression matrix
  target.matrix <- apply(target.matrix, 2, function(x) { (x - mean(x,na.rm = T)) / sd(x,na.rm = T) } )
  input.matrix <- apply(input.matrix, 2, function(x) { (x - mean(x,na.rm =T)) / sd(x,na.rm=T) } )
  input.matrix <- input.matrix[,!is.na(colSums(input.matrix))]
  # setup weight matrix
  num.samples <- dim(target.matrix)[1]
  num.targets <- dim(target.matrix)[2]
  num.inputs <- dim(input.matrix)[2]
  target.names <- colnames(target.matrix)
  input.names <- colnames(input.matrix)
  #print(input.names)
  
  weight.matrix <- matrix(0.0, nrow=num.targets, ncol=num.inputs)
  rownames(weight.matrix) <- target.names
  colnames(weight.matrix) <- input.names
  
  # set mtry
  if (class(K) == "numeric") {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.inputs))
  } else if (K == "all") {
    mtry <- num.inputs-1
  } else {
    stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
  }
  if (trace) {
    cat(paste("Starting RF computations with ", nb.trees,
              " trees/target gene,\nand ", mtry,
              " candidate input genes/tree node\n",
              sep=""))
    flush.console()
  }
  
  # compute importances for every target gene
  names(target.names)<-target.names
  #return(list(target.names,input.matrix,target.matrix))
  numcores = detectCores()
  clst<-makeCluster(numcores)
  clusterExport(cl = clst, varlist = c("importance","randomForest","RSGWM2","num.targets","target.names","input.matrix","target.matrix","trace","mtry","nb.trees","importance.measure",...),envir = environment())
  imList<-parLapplyLB(cl=clst, X=target.names, function(x) RSGWM2(x,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...))
  stopCluster(cl=clst)
  
  #imList<-lapply(target.names,function(x) RSGWM2(x,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...))
  #return(imList)
  for(nm in names(imList))
  {
    tcols<-names(imList[[nm]])
    weight.matrix[nm,tcols] <- imList[[nm]]
  }
  # weight.matrix<-sapply(imList,function(x) x)
  
  mynet <- weight.matrix/num.samples
  mynet <- (mynet-min(mynet,na.rm=TRUE))/(max(mynet,na.rm=TRUE)-min(mynet,na.rm=TRUE))
  return(mynet)
  #    return(list(Weight=weight.matrix,Model=model.matrix,PedictionCorrelations=cor.vec))
}   

RSGWM2<-function(target.gene.name,num.targets,target.names,input.matrix,target.matrix,trace,mtry,nb.trees,importance.measure,...)
{
  target.gene.idx<-which(target.names==target.gene.name)
  if (trace) 
  {
    cat(paste("Computing gene ", target.gene.idx, "/", num.targets, "\n", sep=""))
    flush.console()
  }
  #target.gene.name <- target.names[target.gene.idx]
  # remove target gene from input genes
  #these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
  temp.input.matrix<-input.matrix
  if(target.gene.name %in% colnames(input.matrix))
  {
    #rmind<-which(colnames(input.matrix)==target.gene.name)
    rmind<-grep(target.gene.name,colnames(input.matrix),fixed=T)
    temp.input.matrix<-input.matrix[,-rmind]
    print(paste("Removing:",target.gene.name,"fom Input| Index Number:",rmind," |  New Dimensions:",dim(temp.input.matrix)[1],"X",dim(temp.input.matrix)[2]))
  }
  x <- temp.input.matrix
  y <- target.matrix[,target.gene.name]
  #incSamp<-names(y)[which(!is.na(y))]
  #x<-x[incSamp,]
  #y<-y[incSamp]
  
  rf <- randomForest(x = x, y = y, mtry=mtry, ntree=nb.trees, keep.forest=F, importance=TRUE,...)
  im <- importance(rf)[,importance.measure]

  #im.names <- names(im)
  return(im)
}