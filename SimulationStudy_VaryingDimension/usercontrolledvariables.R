## This file accompanies the manuscript:  Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje "Rotation Invariance in non--Brownian motion Phylogenetic Comparative Methods: A comment on Adams and Collyer (2018)"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## This script contains all the variables with which the user can steer the analysis.


## User controlled variables ===========================================================================
simulsetup_script<-"simulsetups.R" ## file where a list with the modelsetups is located
c_dirpreffix<-"" ## if the results and summaries are (to be) in some other directory than the current one
c_file_suffix<-"" ## suffix of directories and files, can be used for distinguishing between numerous reruns
c_file_suffix_estimation<-"" ## suffix of directories and files, can be used for distinguishing between numerous reruns
c_simuldir<-paste0("./SimulationRuns",c_file_suffix) ## the main directory where the results of the simulations are located 
c_simuldir_estimation<-paste0("./SimulationRuns",c_file_suffix_estimation) ## the main directory where the results of the simulations are located 
c_fileprefix<-"SimulationReestimation_SetupID" ## how the file names with the simulation results begin
c_boxplotdir<-"./BoxPlots" ## directry into which to save the boxplots
c_file_modelcomparison<-paste0("modelcomparisons",c_file_suffix,".txt") ## text file where the model comparisons (by AICc) will be saved
estimation_directory<-paste0(c_simuldir_estimation,"/Estimation",c_file_suffix_estimation,"/")
v_p<-c(4,8,12,16,32) ## dimensions (traits) considered
vN<-c(32,64,128,256,512,1024,2048)
viter<-list("size_32"=1:100,"size_64"=1:100,"size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100)
numcores<-48 ## number of cores used for the analyses
b_should_random_seed_be_saved<-TRUE ## should the random number generator seeds be saved during the analyses rerun
b_use_random_seed_from_manuscript<-FALSE ## should a new analyses be done or from the random seeds in the manuscript (actually whatever ones are in the directories ./SimulationRuns/*_ModelSelection/individual_runs_results_seeds/, where * is the given setup id, in the manuscript these were 01, 02, 03, 04, 05

## list of models under which the data was simulated, in the case of A we have to consider for correct comparison
## two cases diagA positive and NULL, for simulation of course this does not matter, but
## for estimation it does, this list has to be consitent with the lmodelcomparison_setups list found in the simulsetup_script R file
lsimmodellist<-list(
list(id="01",lmodelcf=list(evolmodel="bm_indep"),modelparams=list(),modelnameforplot="bm_indep_4",mvSLOUCHsummary=list()),
list(id="02",lmodelcf=list(evolmodel="bm_indep"),modelparams=list(),modelnameforplot="bm_indep_8",mvSLOUCHsummary=list()),
list(id="03",lmodelcf=list(evolmodel="bm_indep"),modelparams=list(),modelnameforplot="bm_indep_12",mvSLOUCHsummary=list()),
list(id="04",lmodelcf=list(evolmodel="bm_indep"),modelparams=list(),modelnameforplot="bm_indep_16",mvSLOUCHsummary=list()),
list(id="05",lmodelcf=list(evolmodel="bm_indep"),modelparams=list(),modelnameforplot="bm_indep_32",mvSLOUCHsummary=list())
)
vsetupstorun<-c("01"=TRUE,"02"=TRUE,"03"=TRUE,"04"=TRUE,"05"=TRUE) ## which model setups to run in case we would have many

## end of user controlled variables
## =====================================================================================================

# for short testing
# v_p<-c(4,5)#,8,12,16,32) ## dimensions (traits) considered
# vN<-c(5,6)
# viter<-list("size_5"=1:3,"size_6"=1:3)
# lsimmodellist<-list(
# list(id="01",lmodelcf=list(evolmodel="bm_indep"),modelparams=list(),modelnameforplot="bm_indep_4",mvSLOUCHsummary=list()),
# list(id="02",lmodelcf=list(evolmodel="bm_indep"),modelparams=list(),modelnameforplot="bm_indep_8",mvSLOUCHsummary=list())
# )
# vsetupstorun<-c("01"=TRUE,"02"=TRUE) ## which model setups to run in case we would have many

