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

## User controlled variables ===========================================================================
c_dirpreffix<-"" ## if the results and summaries are (to be) in some other directory than the current one
c_simuldir<-"./SimulationRuns" ## the main directory where the results of the simulations are located 
c_fileprefix<-"SimulationReestimation_SetupID" ## how the file names with the simulation results begin
c_file_suffix<-"" ## suffix of directories and files, can be used for distinguishing between numerous reruns, corresponds to runnum in tosource_simulationreestimation.R
c_boxplotdir<-"./BoxPlots" ## directry into which to save the boxplots
c_file_modelcomparison<-paste0("modelcomparisons",c_file_suffix,".txt") ## text file where the model comparisons (by AICc) will be saved
vN<-c(32,64,128,256,512,1024) ## number of tips
v_setups<-c(1,2,3,4,5) ## setup id, the preceding 0 will be added later on in the scripts
## end of user controlled variables
## =====================================================================================================

## correct directory names with suffixes
c_simuldir<-paste0(c_simuldir,c_file_suffix)
c_boxplotdir<-paste0(c_boxplotdir,c_file_suffix)
## =====================================================================================================

## =====================================================================================================
library(mvSLOUCH)
library(TreeSim)
library(ape)
source("summarizeRes.R");
source("simulsetups.R")
## =====================================================================================================


## list of models under which the data was simulated, in the case of A we have to consider for correct comparison
## two cases diagA positive and NULL, for simulation of course this does not matter, but
## for estimation it does
lsimmodellist<-list(
list(id="01",lmodelcf=list(evolmodel="bm"),modelparams=list(),modelnameforplot="BM",mvSLOUCHsummary=list()),
list(id="02",lmodelcf=list(evolmodel="ouch",Atype="Diagonal",Syytype="Diagonal",diagA=NULL),modelparams=list(),modelnameforplot="OUOUs1",mvSLOUCHsummary=list()),
list(id="02",lmodelcf=list(evolmodel="ouch",Atype="Diagonal",Syytype="Diagonal",diagA="Positive"),modelparams=list(),modelnameforplot="OUOUs1P",mvSLOUCHsummary=list()),
list(id="03",lmodelcf=list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA=NULL),modelparams=list(),modelnameforplot="OUOUs2",mvSLOUCHsummary=list()),
list(id="03",lmodelcf=list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA="Positive"),modelparams=list(),modelnameforplot="OUOUs2P",mvSLOUCHsummary=list()),
list(id="04",lmodelcf=list(evolmodel="mvslouch",Atype="Diagonal",Syytype="Diagonal",diagA=NULL),modelparams=list(),modelnameforplot="OUBMs1",mvSLOUCHsummary=list()),
list(id="04",lmodelcf=list(evolmodel="mvslouch",Atype="Diagonal",Syytype="Diagonal",diagA="Positive"),modelparams=list(),modelnameforplot="OUBMs1P",mvSLOUCHsummary=list()),
list(id="05",lmodelcf=list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA=NULL),modelparams=list(),modelnameforplot="OUBMs2",mvSLOUCHsummary=list()),
list(id="05",lmodelcf=list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA="Positive"),modelparams=list(),modelnameforplot="OUBMs2P",mvSLOUCHsummary=list())
)


sink(c_file_modelcomparison)
for (N in vN){for(i in v_setups){
sum_modelcomp<-f_getres(paste0(c_dirpreffix,c_simuldir,"/","0",i,"_ModelSelection/",c_fileprefix,"_0",i,"_N_",N,c_file_suffix,".RData"),lmodelcomparison_setups[[i]]$recoverunder)
cat(paste0("N=",N," setup: 0",i,"\n"))
print(sum_modelcomp$mResSummary)
cat("=====================================================================================================\n")
}}
sink()

source("paramest.R") ## create boxplots
