install.packages("ape")
install.packages("pegas")
install.packages("seqinr")
install.packages("ggplot2")
install.packages("adegenet")
install.packages("hierfstat")

#remove everything
rm(list=ls())


#load libraries
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("hierfstat")

#set working directory
#setwd("C:/Users/UsuarioPC/Desktop/..../")

#where am I working?
#getwd()

#Open STRUCTURE file (or Genepop file)
myFile <- import2genind("Xr-noSNI03-islands.gen")


########################
X <- scaleGen(myFile, NA.method="mean")
pca1<-dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)


myCol <-c("lightblue","darkblue","darkgreen", "black") ## Xantusia all pops
#myCol <-c("lightblue","darkblue","darkgreen") ## Xantusia island pops
#myCol <-c("black","firebrick4","darkorange", "coral2") ## Pseudacris all pops



###GRAPHING PCA####

s.class(pca1$li,pop(myFile),xax=1,yax=2,col=myCol,axesell=FALSE,
        cstar=0,cpoint=1.5,grid=FALSE)



####################################
#####   DAPC by original pops   ####
####################################


dapc2<-dapc(X,pop(myFile))
dapc2 #results dapc
scatter(dapc2, col=myCol) #graph dapc
summary(dapc2) #result summary dapc
contrib<-loadingplot(dapc2$var.contr,axis=1,thres=.07,lab.jitter=1) #contributions of each SNP to linear discriminants




