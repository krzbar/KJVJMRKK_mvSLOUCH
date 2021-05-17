## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## This R script demonstrates that the likelihood calculations of mvSLOUCH are rotation invariant.

library(mvSLOUCH)
library(mvtnorm)
#set.seed(12345)
RNGversion(as.character(getRversion()))
set.seed(1234, kind = "Mersenne-Twister", normal.kind = "Inversion")

#numtips<-30
numtips<-2
numtraits<-2
phyltree<-ape::rtree(numtips)
OUOUdata <- rmvnorm(n = 2, mean = rep(0, 2.5), sigma = matrix(c(1.4,-1,-1,1.4), ncol=2) * 1.5)
rownames(OUOUdata) = rev(phyltree$tip.label)
colnames(OUOUdata) = paste0("X", 1:2)
      
phyltree<-phyltree_paths(phyltree) ## a technical step to speed up calculations, this enhaces the phylogeny with information that will be later used by mvSLOUCH
#OUOUparameters<-list(vY0=matrix(c(0,0,0),nrow=3,ncol=1),A=rbind(c(9,5,0),c(5,5,0),c(1,2,1)),mPsi=matrix(c(0,0,0),nrow=3,ncol=1),Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))
##OUOUparameters<-list(vY0=matrix(c(0,0,0),nrow=3,ncol=1),A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=matrix(c(0,0,0),nrow=3,ncol=1),Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))
#OUOUparameters<-list(vY0=matrix(c(1.2,-1.2),nrow=numtraits,ncol=1),A=rbind(c(3,0.8),c(0.8,3)),mPsi=matrix(c(0.5,0.6),nrow=numtraits,ncol=1),Syy=rbind(c(1,-0.98),c(0,1)))
OUOUparameters<-list(vY0=matrix(c(-1.2,2),nrow=numtraits,ncol=1),A=rbind(c(3,0.8),c(0.8,2)),mPsi=matrix(c(1,1.3),nrow=numtraits,ncol=1),Syy=rbind(c(1.4,-0.6),c(0,1.4)))

##OUOUdata<-mvSLOUCH::simulOUCHProcPhylTree(phyltree,OUOUparameters,NULL,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

OUOU.summary<-mvSLOUCH::SummarizeOUCH(phyltree,OUOUdata,OUOUparameters,NULL,t=c(1),dof=10)

print(OUOU.summary[[1]]$LogLik)

## do the non-phylogenetic "PCA"
#pcares<-stats::prcomp(OUOUdata, center = FALSE, scale. = FALSE)
#U<-pcares$rotation ## rotation matrix to obtain the "PCs", i.e. each observation row of data is multiplied by U
M<-matrix(rnorm(numtraits*numtraits),numtraits,numtraits);
M<-M%*%t(M)
U<-eigen(M)$vectors
U<-M
##OUOUdata_rot<-pcares$x
OUOUdata_rot<-matrix(NA,ncol=numtraits,nrow=numtips)
rownames(OUOUdata_rot)<-rownames(OUOUdata)
colnames(OUOUdata_rot)<-colnames(OUOUdata)
for (i in 1:numtips){
    OUOUdata_rot[i,]<-t(U)%*%OUOUdata[i,]
}

## transform estimated parameters so that they correspond to the rotated data
rot_par<-OUOUparameters
rot_par$A<-(t(U)%*%rot_par$A%*%solve(t(U)))
rot_par$Syy<-(t(U)%*%rot_par$Syy)
rot_par$vY0<-(t(U)%*%rot_par$vY0)
rot_par$mPsi<-(t(U)%*%rot_par$mPsi)



OUOU_rot_summary<-mvSLOUCH::SummarizeOUCH(phyltree,OUOUdata_rot,rot_par,NULL,t=c(1),dof=10)
##OUOU_rot_summary<-mvSLOUCH::SummarizeOUCH(phyltree,OUOUdata_rot,OUOUparameters,NULL,t=c(1),dof=10)
print(OUOU_rot_summary[[1]]$LogLik) ## likelihood should be equal to that of unrotated data
print(OUOU_rot_summary[[1]]$LogLik+numtips*log(abs(det(U)))) ## need to correct for by Jacobian!
