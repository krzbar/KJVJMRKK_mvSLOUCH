## This file accompanies the manuscript: Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje "Analytical advances alleviate model misspecification in non--Brownian multivariate comparative methods"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

f_make_string_from_setup<-function(lsetup){
    simulated_model_name<-"bm"
    if (lsetup$evolmodel=="bm_indep"){simulated_model_name<-"bm_indep"}
    if (!(lsetup$evolmodel%in% c("bm","bm_indep"))){
	## Here column names can become ambigous if diagA is involved    
	simulated_model_name<-paste(lsetup$evolmodel,lsetup$Atype,lsetup$Syytype,sep="_")
	if (!is.null(lsetup$diagA)){
	    simulated_model_name<-paste(simulated_model_name,lsetup$diagA,sep="_")
	}
	if (!is.null(lsetup$estimateBmethod)){
	    simulated_model_name<-paste(simulated_model_name,lsetup$estimateBmethod,sep="_")
	}
	if (!is.null(lsetup$parameter_signs)){
	    for (x in names(lsetup$parameter_signs)){
		simulated_model_name<-paste(simulated_model_name,x,sep="_")
		simulated_model_name<-paste(simulated_model_name,paste(lsetup$parameter_signs[[x]],collapse="_"),sep="_")
	    }
	}
    }
    simulated_model_name
}

extract_estimation_summary<-function(simulsetup,model_setups,mvsl_estim_res){
#save(simulsetup,model_setups,mvsl_estim_res,file="tmpsetups.RData")
    vAICcDiff<-rep(Inf,length(model_setups))
    vBestAICcIndex<-rep(NA,length(model_setups))
    names(vAICcDiff)<-sapply(model_setups,f_make_string_from_setup,simplify=TRUE)
    names(vBestAICcIndex)<-names(vAICcDiff)
    for (i in 1:length(mvsl_estim_res$testedModels)){
	model_type<-f_make_string_from_setup(mvsl_estim_res$testedModels[[i]]$model)
	aic.c_value<-Inf
	if (is.element("aic.c",names(mvsl_estim_res$testedModels[[i]]))){
	    aic.c_value<-mvsl_estim_res$testedModels[[i]]$aic.c
	}
	if (aic.c_value < vAICcDiff[model_type]){
    	    vAICcDiff[model_type]<-aic.c_value
    	    vBestAICcIndex[model_type]<-i
	}
    }    
    simulated_model_name<-f_make_string_from_setup(simulsetup)
    vAICcTrueChosen<-rep(NA,length(vAICcDiff))
    names(vAICcTrueChosen)<-names(vAICcDiff)

    simulated_reestim_results<-NA
    v_simulated_under_models<-startsWith(names(vBestAICcIndex),simulated_model_name)
    if (!all(is.na(vBestAICcIndex[v_simulated_under_models]))){
	#vAICc<-sapply(vBestAICcIndex[v_simulated_under_models],function(i,ltestmodels){ltestmodels[[i]]$aic.c},ltestmodels=mvsl_estim_res$testedModels,simplify=TRUE)
	simunderbestmodel_index<-NA
	bestaicc_index<-NA
	best_index<-NA
	best_aicc<-NA
	for (j in 1:length(vBestAICcIndex)){
	    if (v_simulated_under_models[j]){
		j_index<-vBestAICcIndex[j]
		if (is.na(best_index)){bestaicc_index<-j;best_index<-j_index;best_aicc<-mvsl_estim_res$testedModels[[j_index]]$aic.c}
		else{
		    if (!is.na(j_index)){
			if (best_aicc<mvsl_estim_res$testedModels[[j_index]]$aic.c){
			    bestaicc_index<-j
			    best_index<-j_index
			    best_aicc<-mvsl_estim_res$testedModels[[j_index]]$aic.c	
			}
		    }
		}
	    }
	}
	simunderbestmodel_index<-best_index[[1]]
	simulated_reestim_results<-mvsl_estim_res$testedModels[[simunderbestmodel_index]]$result
        vAICcDiff<-vAICcDiff-vAICcDiff[bestaicc_index]
	vAICcTrueChosen<-(vAICcDiff>0)
	names(vAICcTrueChosen)<-names(vAICcDiff)
	vAICcTrueChosen[bestaicc_index]<-NA    
    }else{
	stop("Estimation failed for the model under which simulation took place")
    }
    lres<-list(
	AICcDiff=vAICcDiff,TrueChosen=vAICcTrueChosen,estimation_under_true=simulated_reestim_results,BestModel=mvsl_estim_res$BestModel$BestModel,TruePoint_model=simulsetup$modelparams
    )    
    lres
}
