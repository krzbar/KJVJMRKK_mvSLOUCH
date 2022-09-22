## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Model Selection Performance in Phylogenetic Comparative Methods under multivariate Ornstein–Uhlenbeck Models of Trait Evolution"

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
ferula_estimation_directory<-paste0("./SimulationRuns_FerulaTree",runnum,"/Ferula_Estimation",runnum,"/")
main_directory<-paste0("./SimulationRuns_FerulaEstimates",runnum,"/")
dir.create(main_directory, showWarnings = FALSE)
vsetupstorun<-c("01"=TRUE,"02"=TRUE,"03"=TRUE,"04"=TRUE,"05"=TRUE,"06"=TRUE,"07"=TRUE,"08"=TRUE,"09"=TRUE,"10"=TRUE)
if (!exists("cl")){cl <- makeCluster(getOption("cl.cores", numcores),outfile="")}

ferula_ME_data_file<-"Data_ME.csv"
# =====================================================

### Read files with Ferula data to get measurement error model ###
feruladata_ME_all <-read.csv(ferula_ME_data_file,header=TRUE,row.names=1)
v_vars<-colnames(feruladata_ME_all)[6:10]
mFerula_Merror<-feruladata_ME_all[,v_vars,drop=FALSE]
save(mFerula_Merror,file=paste0(main_directory,"Ferula_Merror.RData"))
## =====================================================================================================================
source("simulsetups_FerulaEstimates.R")
source("do_simulationstudy_KnownTreeMerrorPossible.R")

