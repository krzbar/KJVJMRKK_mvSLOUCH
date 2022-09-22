## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Model Selection Performance in Phylogenetic Comparative Methods under multivariate Ornstein–Uhlenbeck Models of Trait Evolution"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## This R script performs the analysis of Feeding styles and oral morphology in ungulates under the 
## B&M regime setup.

library(geiger)
library(PCMBaseCpp)
library(mvSLOUCH)

RNGversion("4.1.1")
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

# Modify MP reconstruction based on SCM
SCM.reg<-regimesFitch$branch_regimes
SCM.reg[c(4:6,49,50,66,67,75,85,96:98,103:106,111,117,121,138,154,155,159)]<-"M"
SCM.reg[30]<-"G"
# Clumping M & B
SCM.reg[SCM.reg=="M"]<-"B"
# Confirm
SCM.col<-SCM.reg
SCM.col[SCM.col=="B"]<-"green"
SCM.col[SCM.col=="G"]<-"blue"
plot(PolyTree,cex = 0.3,edge.color = SCM.col)
add.scale.bar(length=0.2,cex=0.5)


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































