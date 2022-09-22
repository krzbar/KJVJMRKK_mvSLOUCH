## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Model Selection Performance in Phylogenetic Comparative Methods under multivariate Ornstein–Uhlenbeck Models of Trait Evolution"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(PCMBaseCpp)
library(mvSLOUCH)
library(ape)

edge_cutoff<-0.02

fprep_edges<-function(phyltree,edge_cutoff,new_height=1){
    phyltree$edge.length<-new_height*(phyltree$edge.length/max(ape::node.depth.edgelength(phyltree)))
    phyltree$edge.length[which(phyltree$edge.length<edge_cutoff)]<-edge_cutoff
    if (!ape::is.ultrametric(phyltree)){
        vnode_heights_all<-node.depth.edgelength(phyltree)
        curr_height<-max(vnode_heights_all)    

	## code from phyltree_paths.R
	vnodes_edgecol1<-sort(unique(phyltree$edge[,1]))
	vnodes_edgecol2<-sort(unique(phyltree$edge[,2]))
	tip_species_index<-sort(setdiff(vnodes_edgecol2,vnodes_edgecol1)) ## only those nodes that are an ending
    
	
	edge_lengthen_all<-rep(0,length(phyltree$edge.length))
	vedge_id<-which(is.element(phyltree$edge[,2],tip_species_index))
	vedge_id<-vedge_id[order(phyltree$edge[vedge_id,2])]
	edge_lengthen_all[vedge_id]<-curr_height-vnode_heights_all[tip_species_index]
	
	phyltree$edge.length<-phyltree$edge.length+edge_lengthen_all
	phyltree$edge.length<-new_height*(phyltree$edge.length/max(ape::node.depth.edgelength(phyltree)))
    }    
    phyltree
}

### Read files ###

ferulatreeape<-read.tree("Ferula_fruits_tree.txt")
ferulatreeape<-fprep_edges(ferulatreeape,edge_cutoff,new_height=1)

feruladata_all <-read.csv("Data_ME.csv",header=TRUE,row.names=1)
feruladata_all <-feruladata_all[ferulatreeape$tip.label,]

v_means<-colnames(feruladata_all)[1:5]
v_vars<-colnames(feruladata_all)[6:10]

feruladata <- data.matrix(feruladata_all[,v_means, drop=FALSE])
M.error<-sapply(1:nrow(feruladata_all),function(i,feruladata_all,v_vars){x<-feruladata_all[i,v_vars];diag(c(x,0), nrow=5, ncol=5)},feruladata_all=feruladata_all,v_vars=v_vars,simplify=FALSE)

seed1 <- 616 
set.seed(seed1)

print(c("Running starting ouchModel..."))

OUOUmodel_start <- ouchModel(ferulatreeape, feruladata, regimes=NULL, Atype="DecomposablePositive", Syytype="Diagonal", diagA=NULL, maxiter=c(25,50))

seed2 <- 2597
set.seed(seed2)

print(c("Running ouchModel..."))

OUOUmodel <- ouchModel(ferulatreeape, feruladata, regimes=NULL, Atype="Any", Syytype="Diagonal", diagA=NULL, maxiter=c(100,1000), M.error=M.error, start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+"))))

print(c("Running parametric.bootstrap..."))

OUOUbootstrap_model7_sDiagonal_var <- parametric.bootstrap(estimated.model=OUOUmodel, phyltree=ferulatreeape, values.to.bootstrap=c("mPsi", "trait.regression", "corr.matrix", "conditional.corr.matrix", "phyl.halflife"), regimes=NULL, M.error=M.error, predictors=NULL , kY=NULL, numboot=1000, Atype="Any", Syytype="Diagonal", diagA=NULL, parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+"))))

# write the full output to the file

capture.output(OUOUbootstrap_model7_sDiagonal_var, file = "full_output_boot_model7_sDiagonal_var.txt")

# CI for trait regression

y <- matrix(,nrow=length(OUOUbootstrap_model7_sDiagonal_var$bootstrapped.parameters$trait.regression), ncol=length(OUOUmodel$FinalFound$ParamSummary$trait.regression[[1]]))

CI <- matrix(,nrow=2, ncol=length(OUOUmodel$FinalFound$ParamSummary$trait.regression[[1]]))
rownames(CI) <- c("BTL", "BTU")

Output_tr <- vector(mode = "list", length = length(OUOUmodel$FinalFound$ParamSummary$trait.regression))

for(i in 1:length(OUOUmodel$FinalFound$ParamSummary$trait.regression)){
       x <- sapply(OUOUbootstrap_model7_sDiagonal_var$bootstrapped.parameters$trait.regression,function(x) x[i])
        
        for(j in 1:length(x)) {y[j,] <- x[[j]][1:length(OUOUmodel$FinalFound$ParamSummary$trait.regression[[1]])]}
         
         for(k in 1:length(OUOUmodel$FinalFound$ParamSummary$trait.regression[[1]])) { 
             BTCI <- quantile(y[,k], c(0.025, 0.975))
             CI[1,k] <- BTCI[1]
             CI[2,k] <- BTCI[2]
              }
          colnames(CI) <- colnames(OUOUmodel$FinalFound$ParamSummary$trait.regression[[i]])
          Output_tr[[i]] <- CI
}

print(Output_tr)
capture.output(Output_tr, file = "CI_trait.regression_model7_sDiagonal_var.txt")

# CI for correlation matrix

Output_cm <- vector(mode="list", length=2)
names(Output_cm) <- c("BTL", "BTU")

BTU.CorrMat <- rep(NA,length(as.vector(OUOUmodel$FinalFound$ParamSummary$corr.matrix)))
BTL.CorrMat<-BTU.CorrMat
for(i in 1:length(as.vector(OUOUmodel$FinalFound$ParamSummary$corr.matrix))){
  BT.CorrMat<-quantile(sapply(OUOUbootstrap_model7_sDiagonal_var$bootstrapped.parameters$corr.matrix,function(x) x[i]),c(0.025,0.975))
  BTL.CorrMat[i] <- BT.CorrMat[1]
  BTU.CorrMat[i] <- BT.CorrMat[2]
}
BTL.CorrMat <- matrix(BTL.CorrMat, nrow =
  nrow(OUOUmodel$FinalFound$ParamSummary$corr.matrix))
BTU.CorrMat <- matrix(BTU.CorrMat, nrow =
  nrow(OUOUmodel$FinalFound$ParamSummary$corr.matrix))
dimnames(BTL.CorrMat) <- dimnames(BTU.CorrMat)<-
  list(row.names(OUOUmodel$FinalFound$ParamSummary$corr.matrix),
      colnames(OUOUmodel$FinalFound$ParamSummary$corr.matrix))

Output_cm[[1]] <- BTL.CorrMat
Output_cm[[2]] <- BTU.CorrMat

print(Output_cm)
capture.output(Output_cm, file = "CI_correlation.matrix_model7_sDiagonal_var.txt")

# CI for halflives

Output_hl <- vector(mode="list", length=2)
names(Output_hl) <- c("BTL", "BTU")

BT_halflives <- vector(mode="list", length=length(OUOUbootstrap_model7_sDiagonal_var$bootstrapped.parameters$phyl.halflife))

for(i in 1:length(OUOUbootstrap_model7_sDiagonal_var$bootstrapped.parameters$phyl.halflife)) {
          BT_halflives[[i]] <- OUOUbootstrap_model7_sDiagonal_var$bootstrapped.parameters$phyl.halflife[[i]]$halflives
          }

BTU.Halflives <- rep(NA,length(as.vector(OUOUmodel$FinalFound$ParamSummary$phyl.halflife$halflives)))
BTL.Halflives<-BTU.Halflives
for(i in 1:length(as.vector(OUOUmodel$FinalFound$ParamSummary$phyl.halflife$halflives))){
  BT.Halflives <-quantile(sapply(BT_halflives,function(x) x[i]),c(0.025,0.975))
  BTL.Halflives[i] <- BT.Halflives[1]
  BTU.Halflives[i] <- BT.Halflives[2]
}
BTL.Halflives <- matrix(BTL.Halflives, nrow =
  nrow(OUOUmodel$FinalFound$ParamSummary$phyl.halflife$halflives))
BTU.Halflives <- matrix(BTU.Halflives, nrow =
  nrow(OUOUmodel$FinalFound$ParamSummary$phyl.halflife$halflives))
dimnames(BTL.Halflives) <- dimnames(BTU.Halflives)<-
  list(row.names(OUOUmodel$FinalFound$ParamSummary$phyl.halflife$halflives),
      colnames(OUOUmodel$FinalFound$ParamSummary$phyl.halflife$halflives))

Output_hl[[1]] <- BTL.Halflives
Output_hl[[2]] <- BTU.Halflives

print(Output_hl)
capture.output(Output_hl, file = "CI_halflives_model7_sDiagonal_var.txt")
