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

ferulatreeape<-read.tree("Ferula_fruits_tree.txt")
ferulatreeape<-fprep_edges(ferulatreeape,edge_cutoff,new_height=1)

feruladata<-read.csv("Data_no_ME.csv",header=TRUE,row.names=1)
feruladata<-feruladata[ferulatreeape$tip.label,]

feruladata <- data.matrix(feruladata)

seed1 <- 551 
set.seed(seed1)

print(c("Running starting ouchModel..."))

OUOUmodel_start <- ouchModel(ferulatreeape, feruladata, regimes=NULL, Atype="DecomposablePositive", Syytype="Diagonal", diagA=NULL, maxiter=c(50,100))

seed2 <- 2532
set.seed(seed2)

print(c("Running ouchModel..."))

OUOUmodel <- ouchModel(ferulatreeape, feruladata, regimes=NULL, Atype="Any", Syytype="UpperTri", diagA=NULL, maxiter=c(100,1000), start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+"))))

print(c("Running parametric.bootstrap..."))

OUOUbootstrap_model4_sUpperTri_novar <-parametric.bootstrap(estimated.model=OUOUmodel, phyltree=ferulatreeape, values.to.bootstrap=c("trait.regression", "limiting.trait.regression", "corr.matrix", "conditional.corr.matrix", "phyl.halflife"),regimes=NULL, M.error=NULL, predictors=NULL , kY=NULL, numboot=1000, Atype="Any", Syytype="UpperTri", diagA=NULL, parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+"))))

capture.output(OUOUbootstrap_model4_sUpperTri_novar,file = "output_boot_model7_sUpperTri_novar.txt")
