library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("hierfstat")


myFile <- read.structure("Xr.stru", onerowperind=FALSE, n.ind=101, n.loc=1160, col.lab=1)
##QUESTIONS:
### Which column contains the population factor ('0' if absent)? 
###answer:2
###Which other optional columns should be read (press 'return' when done)? 
###Which row contains the marker names ('0' if absent)? 
###Answer:1

myFile

summary(myFile)
pop(myFile)
names(myFile)

########################################################################
barplot(myFile$Hexp-myFile$Hobs, main="Heterozygosity: expected-observed",
			ylab="Hexp-Hobs")
########################################################################

####convert file from genind to genpop object##
####myFile2<-genind2genpop(myFile)
####myFile2
####summary(myFile2)
####names(myFile2)
####convert file from genind to genpop object##


########################

X <- scaleGen(myFile, NA.method="mean")
class(X)
dim(X)
X[1:5,1:5]
pca1<-dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
s.class(pca1$li,pop(myFile))
add.scatter.eig(pca1$eig[1:20], 3,1,2)
####to plot PCs 1 and 3 
s.class(pca1$li,pop(myFile),xax=2,yax=3,sub="PCA 2-3",csub=2)
Fst(as.loci(myFile))
###to plot with funky colors
s.class(pca1$li,pop(myFile),xax=1,yax=2,col=funky(15),axesell=FALSE,
		cstar=0,cpoint=3,grid=FALSE)
###############
myCol <-c("orange2","darkgreen","black", "darkblue", "red", "gray")
plot(pca1$li, col=myCol, cex=3)


##########################################
###  PAIRWISE Fst WITH BOOTSTRAPPING   ###
##########################################

pairwise.fst(myFile)
replicate(10,pairwise.fst(myFile,pop=sample(pop(myFile))))


t.test(myFile$Hexp,myFile$Hobs,pair=T,var.equal=TRUE,alter="greater")


################################################
###  NEIGHBOR-JOINING TREE COMPARED TO PCA   ###
################################################
library(ape)
tre<-nj(dist(as.matrix(X)))
tre
plot(tre,typ="fan",cex=0.7)


myCol<-colorplot(pca1$li,pca1$li,transp=TRUE,cex=4)
abline(h=0,v=0,col="grey")


plot(tre,typ="fan",show.tip=FALSE)
tiplabels(pch=20,col=myCol,cex=4)



########################################################
###  DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS   ###
########################################################
grp<-find.clusters(X,max.n.clust=40)
names(grp)
grp$size
table(pop(myFile),grp$grp)
dapc1<-dapc(X,grp$grp)
dapc1
scatter(dapc1)
summary(dapc1)
set.seed=(4)
contrib<-loadingplot(dapc1$var.contr,axis=2,thres=.07,lab.jitter=1)

####################
###  COMPOPLOT   ###
####################

compoplot(dapc1,posi="bottomright",lab="",
			ncol=1,xlab="individuals")



###################################################
###  SPATIAL ANALYSIS OF PRINCIPAL COMPONENTS   ###
###################################################

