## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## running this script runs the simulation-reestimation study presented in the article

## User controlled variables ===========================================================================
runnum<-"" ## suffix of directories and files, can be used for distinguishing between numerous reruns
numcores<-10 ## number of cores used for the analyses
b_should_random_seed_be_saved<-TRUE ## should the random number generator seeds be saved during the analyses rerun
b_use_random_seed_from_manuscript<-TRUE ## should a new analyses be done or from the random seeds in the manuscript (actually whatever ones are in the directories ./SimulationRuns/*_ModelSelection/individual_runs_results_seeds/, where * is the given setup id, in the manuscript these were 01, 02, 03, 04, 05
## end of user controlled variables
## =====================================================================================================
source("tosource_for_mvSLOUCH.R")

## Setup for Ferula estimation 
main_directory<-paste0("./SimulationRuns_FerulaTree",runnum,"/")
ferula_estimation_directory<-paste0(main_directory,"Ferula_Estimation",runnum,"/")
dir.create(main_directory, showWarnings = FALSE)
dir.create(ferula_estimation_directory, showWarnings = FALSE)
vsetupstorun<-c("01"=TRUE,"02"=TRUE,"03"=TRUE,"04"=TRUE,"05"=TRUE,"06"=TRUE,"07"=TRUE,"08"=TRUE,"09"=TRUE,"10"=TRUE)
if (!exists("cl")){cl <- makeCluster(getOption("cl.cores", numcores),outfile="")}

source_scripts<-c(
"model4_sDiagonal_startDecomposablePositive.R",
"model4_sDiagonal_startDecomposablePositive_var.R",
"model4_sUpperTri_startDecomposablePositive.R",
"model4_sUpperTri_startDecomposablePositive_var.R",
"model7_sDiagonal_startDecomposablePositive.R",
"model7_sDiagonal_startDecomposablePositive_var.R",
"model7_sUpperTri_startDecomposablePositive.R",
"model7_sUpperTri_startDecomposablePositive_var.R",
"modelbm.R","modelbm_var.R")
ferula_tree_file<-"Ferula_fruits_tree.txt"
ferula_noME_data_file<-"Data_no_ME.csv"
ferula_ME_data_file<-"Data_ME.csv"
# =====================================================



fprep_edges<-function(phyltree,edge_cutoff,new_height=1){
    phyltree$edge.length<-new_height*(phyltree$edge.length/max(ape::node.depth.edgelength(phyltree)))
    phyltree$edge.length[which(phyltree$edge.length<edge_cutoff)]<-edge_cutoff
    if (!ape::is.ultrametric(phyltree)){
        vnode_heights_all<-node.depth.edgelength(phyltree)
        curr_height<-max(vnode_heights_all)    

        ### code from phyltree_paths.R ###

        vnodes_edgecol1<-sort(unique(phyltree$edge[,1]))
        vnodes_edgecol2<-sort(unique(phyltree$edge[,2]))
        tip_species_index<-sort(setdiff(vnodes_edgecol2,vnodes_edgecol1)) ### only those nodes that are an ending
    
        
        edge_lengthen_all<-rep(0,length(phyltree$edge.length))
        vedge_id<-which(is.element(phyltree$edge[,2],tip_species_index))
        vedge_id<-vedge_id[order(phyltree$edge[vedge_id,2])]
        edge_lengthen_all[vedge_id]<-curr_height-vnode_heights_all[tip_species_index]
        
        phyltree$edge.length<-phyltree$edge.length+edge_lengthen_all
        phyltree$edge.length<-new_height*(phyltree$edge.length/max(ape::node.depth.edgelength(phyltree)))
    }    
    phyltree
}

edge_cutoff<-0.02

### Read files with Ferula data ###
ferulatreeape<-ape::read.tree(ferula_tree_file)
ferulatreeape<-fprep_edges(ferulatreeape,edge_cutoff,new_height=1)

feruladata_noME<-read.csv(ferula_noME_data_file,header=TRUE,row.names=1)
feruladata_noME<-feruladata_noME[ferulatreeape$tip.label,]
feruladata_noME <- data.matrix(feruladata_noME)

feruladata_ME_all <-read.csv(ferula_ME_data_file,header=TRUE,row.names=1)
feruladata_ME_all <-feruladata_ME_all[ferulatreeape$tip.label,]
v_means<-colnames(feruladata_ME_all)[1:5]
v_vars<-colnames(feruladata_ME_all)[6:10]
feruladata_ME <- data.matrix(feruladata_ME_all[,v_means, drop=FALSE])
M.error<-sapply(1:nrow(feruladata_ME_all),function(i,feruladata_ME_all,v_vars){x<-feruladata_ME_all[i,v_vars];diag(c(x,0), nrow=5, ncol=5)},feruladata_ME_all=feruladata_ME_all,v_vars=v_vars,simplify=FALSE)
save(ferulatreeape,feruladata_ME,feruladata_noME,M.error,file="Ferula_Data.RData")
## =====================================================================================================================


lFerulaEstimPars<-parSapply(cl,1:10,function(i,ferulatreeape,feruladata_ME,feruladata_noME,M.error,source_scripts,main_directory,ferula_estimation_directory){
    source("tosource_for_mvSLOUCH.R")
    model_script<-source_scripts[i]
    start_time<-proc.time()
    source(model_script)
    end_time<-proc.time()
    run_time<-end_time-start_time
    file_to_save<-paste0(ferula_estimation_directory,model_script,"Data")
    if (i<9){
	save(OUOUmodel,ferulatreeape,feruladata_ME,feruladata_noME,M.error,model_script,run_time=run_time,file=file_to_save) 
	res<-OUOUmodel$FinalFound$ParamsInModel
    }else{
	save(BMmodel,ferulatreeape,feruladata_ME,feruladata_noME,M.error,model_script,run_time=run_time,file=file_to_save) 
	res<-BMmodel$ParamsInModel
    }
    res
},ferulatreeape=ferulatreeape,feruladata_ME=feruladata_ME,feruladata_noME=feruladata_noME,M.error=M.error,source_scripts=source_scripts,main_directory=main_directory,ferula_estimation_directory=ferula_estimation_directory,simplify=FALSE)
save(lFerulaEstimPars,ferulatreeape,feruladata_ME,feruladata_noME,M.error,source_scripts,file=paste0(ferula_estimation_directory,"/BestModels_FerulaEstimation.RData"))
source("simulsetups_FerulaTree.R")
source("do_simulationstudy_KnownTreeMerrorPossible.R")

