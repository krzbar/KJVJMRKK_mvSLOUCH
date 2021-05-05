## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## This R script performs the analysis of Feeding styles and oral morphology in ungulates under the 
## Full regime setup.

library(geiger)
library(PCMBaseCpp)
library(mvSLOUCH)
library(phytools)

RNGversion("4.0.2")
set.seed(5, kind = "Mersenne-Twister", normal.kind = "Inversion")


########################
### DATA PREPARATION ###
########################

# Data
dat<-read.csv("Data.csv",header = T,row.names = 1)

# Tree
Tree<-read.tree(file = "Tree.tre")
Names<-name.check(Tree,dat)
tip<-Names$tree_not_data
# Prune tree
PrunedTree<-drop.tip(Tree,tip)
name.check(PrunedTree,dat)

# Scaling the tree
tree_height<-phyltree_paths(PrunedTree)$tree_height
ScaledTree<-PrunedTree
ScaledTree$edge.length<-ScaledTree$edge.length/tree_height
phyltree_paths(ScaledTree)$tree_height

# Making polytomies from branches with 0 length
PolyTree<-di2multi(ScaledTree)



# Regime
regimes<-dat$FS[order(match(row.names(dat),PolyTree$tip.label))]
regimesFitch<-fitch.mvsl(PolyTree,regimes,root="B",deltran=TRUE)

# Parsimony
reg.col<-regimesFitch$branch_regimes
reg.col[reg.col=="B"]<-"green"
reg.col[reg.col=="M"]<-"orange"
reg.col[reg.col=="G"]<-"blue"
plot(PolyTree,cex = 0.3, edge.color = reg.col)


# Stochastic character mapping (MP-informed root)
Reg<-as.character(dat$FS)
names(Reg)<-row.names(dat)
name.check(PolyTree,Reg)

# Simulate single stochastic character map using empirical Bayes method
aic<-function(logL,k) 2*k-2*logL
aicc<-function(aic,k,N) aic+((2*k*(k+1))/(N-k-1))

# 1. Equal-rates
mtreeER <- make.simmap(PolyTree, Reg, Q="empirical", model = "ER", pi = c(0.5,0.25,0.25))
aicER<-aic(mtreeER$logL,1)
aiccER<-aicc(aicER,1,length(PolyTree$tip.label))

# 2. Symmetrical
mtreeSYM <- make.simmap(PolyTree, Reg, Q="empirical", model = "SYM", pi = c(0.5,0.25,0.25))
aicSYM<-aic(mtreeSYM$logL,sum(upper.tri(mtreeSYM$Q)))
aiccSYM<-aicc(aicSYM,sum(upper.tri(mtreeSYM$Q)),length(PolyTree$tip.label))

# 3. All-rates-different
mtreeARD <- make.simmap(PolyTree, Reg, Q="empirical", model = "ARD", pi = c(0.5,0.25,0.25))
aicARD<-aic(mtreeARD$logL,sum(upper.tri(mtreeARD$Q))*2)
aiccARD<-aicc(aicARD,sum(upper.tri(mtreeARD$Q))*2,length(PolyTree$tip.label))

# Likelihood ratio test (ER vs SYM)
LR1<-2*(mtreeSYM$logL-mtreeER$logL)                                        # Likelihood ratio
p_value1=1-pchisq(LR1,df=attr(aiccSYM,"df")-attr(aiccER,"df"))             # Significance
# Likelihood ratio test (ARD vs SYM)
LR2<-2*(mtreeARD$logL-mtreeSYM$logL)                                       # Likelihood ratio
p_value2=1-pchisq(LR2,df=attr(aiccARD,"df")-attr(aiccSYM,"df"))            # Significance


# SYM map
cols <- setNames(c("green", "blue", "orange"), sort(unique(Reg)))
plotSimmap(mtreeSYM, cols, pts = FALSE, lwd = 3, fsize=0.4, ftype = "i")
add.scale.bar(length=0.2,cex=0.5)
add.simmap.legend(colors = cols, vertical = FALSE, prompt = FALSE,
                  x = 0.02, y = 25.5)

# Description of SCM
describe.simmap(mtreeSYM)


# Multiple mappings under SYM (preferred model)
mtrees <- make.simmap(PolyTree, Reg, Q="empirical", model = "SYM", pi = c(0.5,0.25,0.25), nsim = 500)
XX <- describe.simmap(mtrees)

# Posterior probabilities
plotSimmap(mtrees[[1]], cols, lwd = 3, pts = F, setEnv = T, fsize=0.4, 
           ftype = "i")
nodelabels(pie = XX$ace, piecol = cols, cex = 0.6)
add.scale.bar(length=0.2,cex=0.5)
add.simmap.legend(colors = cols, vertical = FALSE, prompt = FALSE, 
                  x = 0.02, y = 25.5)

XX$ace



# Modify MP object based on SCM
SCM.reg<-regimesFitch$branch_regimes
SCM.reg[c(4:6,49,50,66,67,75,85,96:98,103:106,111,117,121,138,159)]<-"M"
# Confirm
SCM.col<-SCM.reg
SCM.col[SCM.col=="B"]<-"green"
SCM.col[SCM.col=="M"]<-"orange"
SCM.col[SCM.col=="G"]<-"blue"
plot(PolyTree,cex = 0.3,edge.color = SCM.col)
add.scale.bar(length=0.2,cex=0.5)
add.simmap.legend(colors = cols, vertical = T, prompt = FALSE, 
                  x = 0.02, y = 25.5, fsize = 0.5)

# Circular mapping
plot(PolyTree,type="fan",cex=0.7,label.offset=0.01,edge.width=3,
     edge.color=SCM.col)
add.scale.bar(length=0.2,cex=0.7)
add.simmap.legend(colors = cols, vertical = T, prompt = FALSE, 
                  x = -1.3, y = 1.2, fsize = 0.7)



################
### ANALYSIS ###
################


model_setups<-"basic"

# Preparing tree
mvStree<-phyltree_paths(PolyTree)

# Preparing data
mvData<-data.matrix(dat[,c("HM3","MZW","WM3")])
mvData<-log(mvData)


# OUm
estimResults<-estimate.evolutionary.model(mvStree,mvData,regimes = SCM.reg,
                                          root.regime = "B",repeats=5,
                                          model.setups=model_setups,
                                          predictors=c(3),kY=2,doPrint=TRUE,
                                          pESS=NULL,maxiter=c(10,50,100))

# Sort models according to AICc
OUmAIC<-rep(NA,length(estimResults$testedModels))
for (i in 1:length(OUmAIC)) {
  OUmAIC[i]<-estimResults$testedModels[[i]]$aic.c
}
sort(OUmAIC)

# Save output to file
capture.output(estimResults,file = "OUm.txt")


# OU1
OU1Results<-estimate.evolutionary.model(mvStree,mvData,regimes = NULL,
                                        repeats=5,model.setups=model_setups,
                                        predictors=c(3),kY=2,doPrint=TRUE,
                                        pESS=NULL,maxiter=c(10,50,100))

# Sort models according to AICc
OU1AIC<-rep(NA,length(OU1Results$testedModels))
for (i in 1:length(OU1AIC)) {
  OU1AIC[i]<-OU1Results$testedModels[[i]]$aic.c
}
sort(OU1AIC)

# Saving output to file
capture.output(OU1Results,file = "OU1.txt")




##########
### BT ###
##########

OUOUfinal<-estimResults$testedModels[[32]]$result

BT<-parametric.bootstrap(estimated.model=OUOUfinal,phyltree=mvStree,
                         values.to.bootstrap=c("evolutionary.regression","corr.matrix",
                                               "conditional.corr.matrix"),
                         regimes = SCM.reg,root.regime = "B",
                         predictors=c(3),kY=2,numboot=1000,maxiter=c(10,50,100),
                         Atype=estimResults$BestModel$model$Atype,
                         Syytype=estimResults$BestModel$model$Syytype,
                         diagA=estimResults$BestModel$model$diagA)

# Save output to file
capture.output(BT,file = "BT.txt")


# Conditional correlation
EmpU.CondCorr<-rep(NA,length(as.vector(OUOUfinal$FinalFound$ParamSummary$conditional.corr.matrix)))
EmpL.CondCorr<-EmpU.CondCorr
for(i in 1:length(as.vector(OUOUfinal$FinalFound$ParamSummary$conditional.corr.matrix))){
  Emp.CondCorr<-quantile(sapply(BT$bootstrapped.parameters$conditional.corr.matrix,function(x) x[i]),c(0.025,0.975))
  EmpL.CondCorr[i]<-Emp.CondCorr[1]
  EmpU.CondCorr[i]<-Emp.CondCorr[2]
}
EmpL.CondCorr<-matrix(EmpL.CondCorr,nrow = nrow(OUOUfinal$FinalFound$ParamSummary$conditional.corr.matrix))
EmpU.CondCorr<-matrix(EmpU.CondCorr,nrow = nrow(OUOUfinal$FinalFound$ParamSummary$conditional.corr.matrix))
dimnames(EmpL.CondCorr)<-dimnames(EmpU.CondCorr)<-list(row.names(OUOUfinal$FinalFound$ParamSummary$conditional.corr.matrix),colnames(OUOUfinal$FinalFound$ParamSummary$conditional.corr.matrix))


# Correlation matrix
EmpU.CorrMat<-rep(NA,length(as.vector(OUOUfinal$FinalFound$ParamSummary$corr.matrix)))
EmpL.CorrMat<-EmpU.CorrMat
for(i in 1:length(as.vector(OUOUfinal$FinalFound$ParamSummary$corr.matrix))){
  Emp.CorrMat<-quantile(sapply(BT$bootstrapped.parameters$corr.matrix,function(x) x[i]),c(0.025,0.975))
  EmpL.CorrMat[i]<-Emp.CorrMat[1]
  EmpU.CorrMat[i]<-Emp.CorrMat[2]
}
EmpL.CorrMat<-matrix(EmpL.CorrMat,nrow = nrow(OUOUfinal$FinalFound$ParamSummary$corr.matrix))
EmpU.CorrMat<-matrix(EmpU.CorrMat,nrow = nrow(OUOUfinal$FinalFound$ParamSummary$corr.matrix))
dimnames(EmpL.CorrMat)<-dimnames(EmpU.CorrMat)<-list(row.names(OUOUfinal$FinalFound$ParamSummary$corr.matrix),colnames(OUOUfinal$FinalFound$ParamSummary$corr.matrix))


# Evolutionary regression
EmpU.EvoReg<-OUOUfinal$FinalFound$ParamSummary$evolutionary.regression
EmpU.EvoReg[]<-0L
EmpL.EvoReg<-EmpU.EvoReg
for(i in 1:nrow(OUOUfinal$FinalFound$ParamSummary$evolutionary.regression)){
  Emp.EvoReg<-quantile(sapply(BT$bootstrapped.parameters$evolutionary.regression,function(x) x[i]),c(0.025,0.975))
  EmpL.EvoReg[i,]<-Emp.EvoReg[1]
  EmpU.EvoReg[i,]<-Emp.EvoReg[2]
}












































