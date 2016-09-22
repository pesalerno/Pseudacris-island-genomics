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
s.class(pca1$li,pop(myFile), col=myCol)
add.scatter.eig(pca1$eig[1:20], 3,1,2)
####to plot PCs 1 and 3 
s.class(pca1$li,pop(myFile),xax=2,yax=3,sub="PCA 2-3",csub=2)

###to plot with funky colors
s.class(pca1$li,pop(myFile),xax=1,yax=2,col=myCol,axesell=FALSE,
		cstar=0,cpoint=3,grid=FALSE)
###############
myCol <-c("orange2","darkorange","darkgreen", "green2", "darkred", "red2")
plot(pca1$li, col=myCol, cex=3)



##########################################
###         Per locus Fst values       ###
##########################################

Fst(as.loci(myFile))


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
table(pop(),grp$grp)

###dapc by cluster
dapc1<-dapc(X,grp$grp)
dapc1


###dapc by original pop
dapc2<-dapc(X,myFile$pop)
dapc2


scatter(dapc2)
summary(dapc2)

set.seed=(4)
contrib<-loadingplot(dapc1$var.contr,axis=2,thres=.07,lab.jitter=1)

####################
###  COMPOPLOT   ###
####################

compoplot(dapc2,posi="bottomright",lab="",
			ncol=1,xlab="individuals")



###################################################
###  SPATIAL ANALYSIS OF PRINCIPAL COMPONENTS   ###
###################################################





###################################################
###          PCADAPT OUTLIER LOCI TEST          ###
###################################################
library(pcadapt)
###input file formats supported are vcf, ped, lfmm, pcadapt
library(qvalue)


path_to_file<-"./path/file.ped"
myData<-read.pcadapt(path_to_file,type="ped")
print(myData)

##1.perform PCA with large number of PCs (>20)
X<-pcadapt(myData,K=20,transpose=TRUE)


plot(X,option="screeplot")

plot(X,option="scores",pop=poplist)

plot(X,option="scores",i=3,j=4,pop=poplist)

X<-pcadapt(myData,K=2)
summary(X)

qval<-qvalue(X$pvalues)$qvalues
alpha<-0.1
outliers<-which(qval<alpha)
print(outliers)

##which PCs are most correlated to the outliers?
snp_pc<-get.pc(X,outliers)
head(snp_pc)


###################################################
###       POPGENREPORT - POPULATION STATS       ###
###################################################
##popgenreport tutorial general code


library(PopGenReport)
require(PopGenReport)



myFile <- read.genetable("./Xr.stru", sep=" ")
head(platy)

platy.gen <- read.genetable(paste(.libPaths()[1], "/PopGenReport/extdata/platypus1c.csv", sep=""), ind=1, pop=2, lat=3, long=4, other.min=5, other.max=6, oneColPerAll=FALSE, sep="/", ploidy=2)
platy.gen #to check the formatted genind data

##genetic data for each individual is stored within the tab slot
myFile@tab
pop(myFile)
str(myFile@other)
myFileout1<-popgenreport(myFile, mk.counts=TRUE, mk.pdf=FALSE, foldername="test1")
myFileout1
summary(myFileout1)
barplot(myFileout1$counts@nallelesbypop)

complete.test1<-popgenreport(myFile, mk.complete=TRUE, mk.Rcode=TRUE, mk.pdf=FALSE)
summary(complete.test1)
complete.test1
complete.test1$fst ##fst and fis per locus and per population pair
complete.test1$allel.rich
complete.test1$counts
complete.test1$allele.dist ##list of private alleles per population and frequencies

barplot(complete.test1$fst)
plot(complete.test1$FSTbyloci, complete.test1$pop)
scatter(complete.test1$fst$FSTbyloci)

###if you install LaTeX then you can do mk.pdf=TRUE, otherwise it will fail
##platy.out1 <- popgenreport(platy.gen, mk.counts=TRUE, mk.pdf=FALSE, foldername='platy')
##summary(platy.out1)
##barplot(platy.out1$count$nallelesbypop)

