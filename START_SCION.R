#Version 4.2
#Oct 2024

#run all lines in this file to start the SCION RShiny App
#make sure to change this first line to your working directory where SCION lives
setwd('/Users/nclark/Desktop/SCION')

#check that all packages are installed, and load them using pacman
if (!require('pacman', character.only=T, quietly=T)) {
  install.packages('pacman')
  library('pacman', character.only=T)
}else{
  library('pacman',character.only=T)
}

p_load(dtwclust)
p_load(plotly)
p_load(randomForest)
p_load(fastICA)
p_load(shiny)
p_load(shinyscreenshot)
p_load(cluster)
p_load(doParallel)

#for R-4.3 compatibility
library(plotly)
library(randomForest)

source("SCION.R")
source('dtw_clustering.R')
source('ica_clustering.R')
source('kmeans_clustering.R')
source('RS.Get.Weight.Matrix.R')
runApp('app.R')