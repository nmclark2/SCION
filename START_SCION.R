#Version 2.0 beta
#July 2020

#run all lines in this file to start the SCION RShiny App

#check that all packages are installed, and load them
for (package in c('dtwclust', 'plotly','randomForest','fastICA','shiny')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }else{
    library(package,character.only=T)
  }
}

source("SCION.R")
source('dtw_clustering.R')
source('ica_clustering.R')
source('RS.Get.Weight.Matrix.R')
runApp('app.R')