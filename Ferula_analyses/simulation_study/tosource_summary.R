## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## Running this script provides the summary of the resuling simulation-reestimation study
## in particular the boxplots and model comparison abilities
load("Ferula_Data.RData")



simulsetup_script<-"simulsetups_FerulaModelData.R"
## User controlled variables ===========================================================================
c_dirpreffix<-"" ## if the results and summaries are (to be) in some other directory than the current one
c_file_suffix<-"" ## suffix of directories and files, can be used for distinguishing between numerous reruns, corresponds to runnum in tosource_simulationreestimation.R
c_file_suffix_estimation<-"" ## suffix of directories and files, can be used for distinguishing between numerous reruns, corresponds to runnum in tosource_simulationreestimation.R
c_simuldir<-paste0("./SimulationRunsJoint",c_file_suffix) ## the main directory where the results of the simulations are located 
c_simuldir_estimation<-paste0("./SimulationRuns_FerulaTree",c_file_suffix_estimation) ## the main directory where the results of the simulations are located 
c_fileprefix<-"SimulationReestimation_SetupID" ## how the file names with the simulation results begin
c_boxplotdir<-"./BoxPlots" ## directry into which to save the boxplots
c_file_modelcomparison<-paste0("modelcomparisons_FerulaTreeEstimation",c_file_suffix,".txt") ## text file where the model comparisons (by AICc) will be saved
vN<-c(78,128,256,512,1024,2048) ## number of tips
v_setups<-c(1,2,3,4,5,6,7,8,9,10) ## setup id, the preceding 0 will be added later on in the scripts
ferula_estimation_directory<-paste0(c_simuldir_estimation,"/Ferula_Estimation",c_file_suffix_estimation,"/")
num_traits<-5 ## number of traits in the simulation
## end of user controlled variables
## =====================================================================================================

## correct directory names with suffixes
#c_simuldir<-paste0(c_simuldir,c_file_suffix)
c_boxplotdir<-paste0(c_boxplotdir,c_file_suffix)
## =====================================================================================================

## =====================================================================================================
library(mvSLOUCH)
library(TreeSim)
library(ape)
load(paste0(ferula_estimation_directory,"BestModels_FerulaEstimation.RData"))
load(paste0(c_simuldir,"/Ferula_Merror.RData"))
source("summarizeRes.R")
source(simulsetup_script)
## =====================================================================================================


## list of models under which the data was simulated, in the case of A we have to consider for correct comparison
## two cases diagA positive and NULL, for simulation of course this does not matter, but
## for estimation it does
lsimmodellist<-list(
list(id="01",lmodelcf=list(evolmodel="ouch",Atype="Any",Syytype="Diagonal",diagA=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),modelparams=list(),modelnameforplot="m4SD_noME",mvSLOUCHsummary=list()),
list(id="02",lmodelcf=list(evolmodel="ouch",Atype="Any",Syytype="Diagonal",diagA=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),modelparams=list(),modelnameforplot="m4SD_ME",mvSLOUCHsummary=list()),
list(id="03",lmodelcf=list(evolmodel="ouch",Atype="Any",Syytype="UpperTri",diagA=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),modelparams=list(),modelnameforplot="m4UT_noME",mvSLOUCHsummary=list()),
list(id="04",lmodelcf=list(evolmodel="ouch",Atype="Any",Syytype="UpperTri",diagA=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),modelparams=list(),modelnameforplot="m4UT_ME",mvSLOUCHsummary=list()),
list(id="05",lmodelcf=list(evolmodel="ouch",Atype="Any",Syytype="Diagonal",diagA=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),modelparams=list(),modelnameforplot="m7SD_noME",mvSLOUCHsummary=list()),
list(id="06",lmodelcf=list(evolmodel="ouch",Atype="Any",Syytype="Diagonal",diagA=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),modelparams=list(),modelnameforplot="m7SD_ME",mvSLOUCHsummary=list()),
list(id="07",lmodelcf=list(evolmodel="ouch",Atype="Any",Syytype="UpperTri",diagA=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),modelparams=list(),modelnameforplot="m7UT_noME",mvSLOUCHsummary=list()),
list(id="08",lmodelcf=list(evolmodel="ouch",Atype="Any",Syytype="UpperTri",diagA=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),modelparams=list(),modelnameforplot="m7UT_ME",mvSLOUCHsummary=list()),
list(id="09",lmodelcf=list(evolmodel="bm"),modelparams=list(),modelnameforplot="bm_noME",mvSLOUCHsummary=list()),
list(id="10",lmodelcf=list(evolmodel="bm"),modelparams=list(),modelnameforplot="bm_ME",mvSLOUCHsummary=list())
)

sink(c_file_modelcomparison)
for (N in vN){for(i in v_setups){
if (i<10){i0<-paste0("0",i)}
else{i0<-i}
sum_modelcomp<-f_getres(paste0(c_dirpreffix,c_simuldir,"/",i0,"_ModelSelection/",c_fileprefix,"_",i0,"_N_",N,c_file_suffix,".RData"),lmodelcomparison_setups[[i]]$recoverunder)
cat(paste0("N=",N," setup: ",i0,"\n"))
print(sum_modelcomp$mResSummary)
cat("=====================================================================================================\n")
}}
sink()



source("paramest.R") ## create boxplots
