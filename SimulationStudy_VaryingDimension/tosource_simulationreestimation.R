## This file accompanies the manuscript:  Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje "Analytical advances alleviate model misspecification in non--Brownian multivariate comparative methods"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## running this script runs the simulation-reestimation study presented in the article

source("tosource_for_mvSLOUCH.R")

## Setup for estimation 
main_directory<-paste0(c_simuldir_estimation,"/")
runnum<-c_file_suffix
dir.create(main_directory, showWarnings = FALSE)
if (!exists("cl")){cl <- makeCluster(getOption("cl.cores", numcores),outfile="")}
# =====================================================

source(simulsetup_script)
source("do_simulationstudy.R")

