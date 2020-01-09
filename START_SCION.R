#run all lines in this file to start the SCION RShiny App

library(shiny)
setwd("D:/OneDrive/Documents/Walley Lab Postdoc/SCION") #change to your working directory that contains SCION!
source("SCION.R")
source('dtw_clustering.R')
source('RS.Get.Weight.Matrix.R')
runApp('app.R')