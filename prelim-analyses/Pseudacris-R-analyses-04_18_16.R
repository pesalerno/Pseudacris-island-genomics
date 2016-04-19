


******************************************************************
******************************************************************
******************************************************************
#RUNNING PACKAGE: adegenet
******************************************************************
******************************************************************
******************************************************************
library("adegenet")
myFile <- read.structure("Xari-p6-7455snps-95ind-f.stru", onerowperind=FALSE, n.ind=95, n.loc=7455, col.lab=1)
##QUESTIONS:
### Which column contains the population factor ('0' if absent)? 
###answer:2
###Which other optional columns should be read (press 'return' when done)? 
###Which row contains the marker names ('0' if absent)? 
###Answer:1


myFile



*********************************************************************
################PCA with my data......

scaling <- scaleGen(myFile, scale=FALSE)
pca1 <- dudi.pca(scaling, center=FALSE, scale=FALSE) ##displays eigenvalues barplot >>>TELL NUMBER OF RETAINED AXES
pca1


##if you want individual labels in graph
s.label(pca1$li)


		
###Plotting with personalized colors		
myCol <-c("orange2","darkgreen","black", "darkblue", "red", "gray")
s.class(pca1$li, fac=pop(myFile),
		col=myCol, axesel=FALSE, cstar=0, cpoint=1.5)



##############################
*********************************************************************
*********************************************************************
##find.clusters transforms the data using PCA, then it runs k-means algorithm with increasing values of K, and computes summary statistics (default BIC). 
grp <- find.clusters(myFile, max.n.clust=40) #the maximum number of clusters here is k=40
##need to pick number of PCs and number of clusters to retain
grp$size
table(pop(myFile), grp$grp)##how well have groups been retrieved by clustering
table.value(table(pop(myFile), grp$grp), col.lab=paste("inf", 1:6), row.lab=paste("ori", 1:4))
###how many clusters are actually in the data?


********************************************
********************************************
#######DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS ###
dapc <- dapc(myFile)
##don't keep that many PCs, since it may result in instability of membership probabilities. Retain as many as possible without sacrificing much information
dapc

summary(dapc)

###graph with personalized colors:
myCol <-c("orange2","darkgreen","darkblue", "red", "gray")
scatter(dapc, posi.da="bottomright", bg="white",cstar=0, col=myCol)



###COMPUTE CONTRIBUTIONS OF ALLELES IN DAPC
contrib <- loadingplot(dapc$var.contr, axis=2, thres=.07, lab.jitter=1)


*********************************************************************
*********************************************************************
####ASSIGNMENT TO GROUPS AND CLUSTERING
class(dapc$posterior)
dim(dapc$posterior)
round(head(dapc$posterior), 3)
summary(dapc)
assignplot(dapc)

summary(dapc$posterior)
summary(dapc$prior)


******************************************************************
******************************************************************
******************************************************************
#RUNNING PACKAGE: hierfstat
******************************************************************
******************************************************************
******************************************************************
library(hierfstat)

###IMPORTANT#####HIERFSTAT SHOULD BE RUN TWICE, ONCE INCLUDING THE ENTIRE MATRIX, AND ONCE MORE EXCLUDING THE MAINLAND. WHAT IS THE EFFECT OF MAINLAND EXCLUSION WHEN LOOKING AT ISLAND GLOBAL FSTs??


###loading the data into R
data <- read.table("Xari-p6-7455snps-95ind-e.stru", header=TRUE, na.string="0")
##if missing data coded as "0" then it needs to add: na.string="0"
###format of example data file: 
#######eight total columns
############1st columns are the hierarchical population levels
############Last columns are the different loci.
head(data)
##to define which part of the data contains the loci, write:
##69 loci-->>columns[4:72]
loci <- data[,c(4:7455)]
levels <- data[,c(2:3)] #this defines the hierarchical levels for the Fst
head(levels)

###in order to estimate hierarchical Fst from these data, type:
varcomp.glob(levels,loci, diploid=TRUE)

###to test significance of genetic differentiation at the different levels controlling for the effects at the other levels:
test.within(loci,test=pops,within=islands,nperm=1000)
##carry out 1000 permutations of individuals between units defined by pops, but keeping then within units defined by islands

test.between(loci,rand.unit=pops,test=islands,nperm=1000)
##carry out permutations of whole units of lev2 agains units defined by lev1 1000 times

###bootstrapping over loci
boot.vc(levels,loci)



