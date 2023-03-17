## This file accompanies the manuscript:  Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje "Analytical advances alleviate model misspecification in non--Brownian multivariate comparative methods"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## Running this script provides the summary of the resuling simulation-reestimation study
## in particular the boxplots and model comparison abilities

c_boxplotdir<-paste0(c_boxplotdir,c_file_suffix)
v_setups<-1:length(v_p)
num_traits<-v_p
source("summarizeRes.R")
source(simulsetup_script)

sink(c_file_modelcomparison)
for (N in vN){for(i in v_setups){
if (i<10){i0<-paste0("0",i)}
else{i0<-i}
bDoSum<-TRUE
if (bDoSum){
tmp_recoverunder<-lmodelcomparison_setups[[i]]$recoverunder
if (is.element("non_mvSLOUCH_model_setups",names(tmp_recoverunder))&&(!is.null(tmp_recoverunder$non_mvSLOUCH_model_setups))&&(length(tmp_recoverunder$non_mvSLOUCH_model_setups)>0)){
    tmp_recoverunder$model_setups<-c(tmp_recoverunder$model_setups,tmp_recoverunder$non_mvSLOUCH_model_setups)
}

sum_modelcomp<-f_getres(paste0(c_dirpreffix,c_simuldir,"/",i0,"_ModelSelection/",c_fileprefix,"_",i0,"_N_",N,c_file_suffix,".RData"),tmp_recoverunder)
cat(paste0("N=",N," setup: ",i0,"\n"))
print(sum_modelcomp$mResSummary)
cat("=====================================================================================================\n")
}}
}
sink()


source("paramest.R") ## create boxplots
