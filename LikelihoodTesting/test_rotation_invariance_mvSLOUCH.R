## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Rotation Invariance in non–Brownian motion Phylogenetic Comparative Methods: A comment on Adams and Collyer (2018)"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## This R script demonstrates that the likelihood calculations of mvSLOUCH are rotation invariant.

library(mvSLOUCH)

## set.seed(12345) ## if the user wishes to test for a particular random seed
numtips<-30 ## number of tips of the phylogenetic tree, increasing will of course increas the running time
phyltree<-ape::rtree(numtips)
      
phyltree<-phyltree_paths(phyltree) ## a technical step to speed up calculations, this enhaces the phylogeny with information that will be later used by mvSLOUCH
OUOUparameters<-list(vY0=matrix(c(0,0,0),nrow=3,ncol=1),A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=matrix(c(0,0,0),nrow=3,ncol=1),Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

OUOUdata<-mvSLOUCH::simulOUCHProcPhylTree(phyltree,OUOUparameters,NULL,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

OUOU.summary<-mvSLOUCH::SummarizeOUCH(phyltree,OUOUdata,OUOUparameters,NULL,t=c(1),dof=10)
print(OUOU.summary[[1]]$LogLik)

## do the non-phylogenetic "PCA"
pcares<-stats::prcomp(OUOUdata, center = FALSE, scale. = FALSE)
U<-pcares$rotation ## rotation matrix to obtain the "PCs", i.e. each observation row of data is multiplied by U
OUOUdata_rot<-pcares$x

## transform estimated parameters so that they correspond to the rotated data
rot_par<-OUOUparameters
rot_par$A<-(t(U)%*%rot_par$A%*%solve(t(U)))
rot_par$Syy<-(t(U)%*%rot_par$Syy)
rot_par$vY0<-(t(U)%*%rot_par$vY0)
rot_par$mPsi<-(t(U)%*%rot_par$mPsi)



OUOU_rot_summary<-mvSLOUCH::SummarizeOUCH(phyltree,OUOUdata_rot,rot_par,NULL,t=c(1),dof=10)

## likelihood should be equal to that of unrotate data
## need to correct for by Jacobian, if U is orthogonal log(det(U))=0
print(OUOU_rot_summary[[1]]$LogLik+numtips*log(abs(det(U))))


## do same thing for BM model
BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))

BM.summary<-mvSLOUCH::SummarizeBM(phyltree,OUOUdata,BMparameters,t=c(1),dof=10)
print(BM.summary[[1]]$LogLik)
BM_rot<-BMparameters
BM_rot$vX0<-t(U)%*%BM_rot$vX0
BM_rot$Sxx<-t(U)%*%BM_rot$Sxx
BM_rot_summary<-mvSLOUCH::SummarizeBM(phyltree,OUOUdata_rot,BM_rot,t=c(1),dof=10)
print(BM_rot_summary[[1]]$LogLik+numtips*log(abs(det(U)))) ## need to correct for by Jacobian, if U is orthogonal log(det(U))=0


## now let us do estimation to see how the back transformation works
## There is no guarantee of A's class to be retained after the PC transformation
OUOU_PCestimate<-mvSLOUCH::ouchModel(phyltree,OUOUdata_rot,Atype="DecomposablePositive")
OUOU_PCestimate_result<-NA
if (is.list(OUOU_PCestimate$MaxLikFound)){
    OUOU_PCestimate_result<-OUOU_PCestimate$MaxLikFound
}else{OUOU_PCestimate_result<-OUOU_PCestimate$FinalFound}
print(OUOU_PCestimate_result$ParamSummary$LogLik)

unrot_par_est<-OUOU_PCestimate_result$ParamsInModel
unrot_par_est$A<-(solve(t(U))%*%unrot_par_est$A%*%(t(U)))
unrot_par_est$Syy<-(solve(t(U))%*%unrot_par_est$Syy)
unrot_par_est$vY0<-(solve(t(U))%*%unrot_par_est$vY0)
unrot_par_est$mPsi<-(solve(t(U))%*%unrot_par_est$mPsi)

OUOU_unrotest_summary<-mvSLOUCH::SummarizeOUCH(phyltree,OUOUdata,unrot_par_est,dof=OUOU_PCestimate_result$ParamSummary$dof)
## likelihood should be equal to that of unrotated data
## need to correct for by Jacobian, if U is orthogonal log(det(U))=0
## minus here as we take solve(t(U)) not t(U)
print(OUOU_unrotest_summary[[1]]$LogLik-numtips*log(abs(det(U)))) 


## now let us do estimation to see how the back transformation works for BM
BM_PCestimate<-mvSLOUCH::BrownianMotionModel(phyltree,OUOUdata_rot)
print(BM_PCestimate$ParamSummary$LogLik)

unrot_par_est<-BM_PCestimate$ParamsInModel
unrot_par_est$Sxx<-(solve(t(U))%*%unrot_par_est$Sxx)
unrot_par_est$vX0<-(solve(t(U))%*%unrot_par_est$vX0)

BM_unrotest_summary<-mvSLOUCH::SummarizeBM(phyltree,OUOUdata,unrot_par_est,dof=BM_PCestimate$ParamSummary$dof)
print(BM_unrotest_summary[[1]]$LogLik) ## likelihood should be equal to that of unrotated data


