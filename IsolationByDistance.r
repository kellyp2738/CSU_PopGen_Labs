# Isolation by distance
# tested on OSX 30 Sept 2015
# tested on Windows 1 Oct 2015
# key assmuption: the ordering of populations in the geo. and genetic distance matrices are preserved

library(StAMPP)

# load and manipulate the geographic distance
linear<-read.csv('~/Dropbox/BZ 577/Week 5-6/FigsData/linear.csv')
linear<-read.csv('C:/Users/AntolinLab/Dropbox/BZ 577/Week 5-6/FigsData/linear.csv')
linear.mat<-as.matrix(linear[1:13, 2:14])
linear.dist<-log(linear.mat[lower.tri(linear.mat)])

# load and manipulate the genetic distance
# data are automatically loaded as GENIND format
pdogdata<-read.genepop('~/Dropbox/BZ 577/Week 5-6/FigsData/Cynomys population study.gen')
pdogdata<-read.genepop('C:/Users/AntolinLab/Dropbox/BZ 577/Week 5-6/FigsData/Cynomys population study.gen')

# StAMPP needs genlight data... 
pdog.gl<-new("genlight", gen=pdogdata@tab, ploidy=2,
             pop=pdogdata@pop,ind.names=seq(1,155,1), parallel=F)

# calculate Fst matrix
Fstats<-stamppFst(pdog.gl, nboots=10000)

# pull out only the Fst values (also in this data structure are Pvalues and the Fst calculated for each bootstrap iteration)
Fst.mat<-Fstats$Fsts #lower triangle only

# extract the lower triangle
Fst.vec<-Fst.mat[lower.tri(Fst.mat)]

# linearize Fst
Fst.norm<-Fst.vec/(1-Fst.vec)

# visualize data
plot(linear.dist, Fst.norm, xlab='log(Geographic Distance)', 
     ylab='Genetic Distance', pch=16)

# conduct the Mantel test (tell R to use the ecodist version)
PdogMantel <- ecodist::mantel(Fst.norm ~ linear.dist, nperm = 1000)