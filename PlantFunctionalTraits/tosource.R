## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(phytools)
library(geiger)
library(PCMBaseCpp)
library(mvSLOUCH)
library(ape)
library(doParallel)
library(parallel)
if (b_make_optima_plot){
    library(tidyr)
    library(ggplot2)
}

### function to store objects if running analyses in parallel
multiResultClass <- function(random_seed=NULL,mvSLOUCH_output=NULL, runtime=NULL)
{
  me <- list(
    random_seed = random_seed,
    mvSLOUCH_output = mvSLOUCH_output,
    runtime = runtime
  )
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}


f_findBestModel<-function(l_mvSLOUCH_output,information_criterion="aic.c"){
    best_inf_crit_val<-Inf
    best_index<-NA
    best_output<-NA
    for (i in 1:length(l_mvSLOUCH_output)){
	if (is.element("mvSLOUCH_output",names(l_mvSLOUCH_output[[i]]))){mvSLOUCH_output<-l_mvSLOUCH_output[[i]]$mvSLOUCH_output}
	else{mvSLOUCH_output<-l_mvSLOUCH_output[[i]]}
	if (is.element("FinalFound",names(mvSLOUCH_output))){
	    if ((!is.list(mvSLOUCH_output$MaxLikFound))&&(mvSLOUCH_output$MaxLikFound[1]=="Same as final found")){
		output_to_check<-mvSLOUCH_output$FinalFound
	    }else{
		if (is.list(mvSLOUCH_output$MaxLikFound)&&(is.element("ParamSummary",names(mvSLOUCH_output$MaxLikFound)))){
		    output_to_check<-mvSLOUCH_output$MaxLikFound
		}else{output_to_check<-NA}		
	    }	
	}else{
	    if (is.element("ParamSummary",names(mvSLOUCH_output))){
		output_to_check<-mvSLOUCH_output
	    }else{output_to_check<-NA}
	}
	if ((is.list(output_to_check))&&(!is.na(output_to_check[1]))){
	    curr_inf_crit_val<-NA
	    if(information_criterion %in% c("R2","R2_phylaverage")){
	    	infcrit_index<-which(names(output_to_check$ParamSummary$RSS)==information_criterion)
		if(length(infcrit_index)==1){
	    	    curr_inf_crit_val<-output_to_check$ParamSummary$RSS[[infcrit_index]] ## this line of code will need to be updated if other versions of R2 are considered
		    curr_inf_crit_val<-curr_inf_crit_val*(-1) ## R2 needs to be maximized
		}
	    }else{
		infcrit_index<-which(names(output_to_check$ParamSummary)==information_criterion)
		if(length(infcrit_index)==1){curr_inf_crit_val<-output_to_check$ParamSummary[[infcrit_index]]}
	    }
	    if(!is.na(curr_inf_crit_val)){
		if (curr_inf_crit_val<best_inf_crit_val){
		    best_inf_crit_val<-curr_inf_crit_val
		    best_index<-i
		    best_output<-output_to_check
		}
	    }
	}
    }
    list(best_output=best_output,best_index=best_index,inf_crit=information_criterion,inf_crit_value=best_inf_crit_val)
}

f_extractBestFromPair<-function(model_list1,model_list2){
    best_model1<-f_findBestModel(model_list1)
    best_model2<-f_findBestModel(model_list2)
    best_model<-NA
    best_model_startingpoint<-NULL
    best_model_value<-Inf
    best_index<-NA
    which_of_pair<-NA
    if ((!is.na(best_model1$inf_crit_value))&&(!is.null(best_model1$inf_crit_value))&&(!is.infinite(best_model1$inf_crit_value))){
	best_model<-best_model1$best_output
	best_model_value<-best_model1$inf_crit_value
	which_of_pair<-1
	best_index<-best_model1$best_index
    }
    if ((!is.na(best_model2$inf_crit_value))&&(!is.null(best_model2$inf_crit_value))&&(!is.infinite(best_model2$inf_crit_value))){
	if (best_model2$inf_crit_value<best_model_value){
    	    best_model<-best_model2$best_output
    	    best_model_value<-best_model2$inf_crit_value
	    which_of_pair<-2
	    best_index<-best_model2$best_index
	}
    }
    if (is.list(best_model)&&(!is.na(best_model[1]))){
	best_model_startingpoint<-list(A=best_model$ParamsInModel$A,Syy=best_model$ParamsInModel$Syy)
    }
    list(best_model_startingpoint=best_model_startingpoint,best_model_value=best_model_value,best_model=best_model,which_of_pair=which_of_pair,best_index=best_index)
}

f_createstartingpoint<-function(mvSLOUCH_output){
    if (is.element("FinalFound",names(mvSLOUCH_output))){
	    if ((!is.list(mvSLOUCH_output$MaxLikFound))&&(mvSLOUCH_output$MaxLikFound[1]=="Same as final found")){
		output_to_extract_from<-mvSLOUCH_output$FinalFound
	    }else{
		if (is.list(mvSLOUCH_output$MaxLikFound)&&(is.element("ParamSummary",names(mvSLOUCH_output$MaxLikFound)))){
		    output_to_extract_from<-mvSLOUCH_output$MaxLikFound
		}else{output_to_extract_from<-NA}		
	    }	
    }
    l_startinpoint<-NULL
    if (is.list(output_to_extract_from)&&(!is.na(output_to_extract_from))){
	l_startinpoint<-list(A=output_to_extract_from$ParamsInModel$A,Syy=output_to_extract_from$ParamsInModel$Syy)
    }
}

