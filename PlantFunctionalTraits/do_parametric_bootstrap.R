## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

load(paste0("./",c_resultsdirectory,"/model1/OUOU.model1.BestModel.RData"))
load(paste0("./",c_resultsdirectory,"/model2/OUOU.model2.BestModel.RData"))
load(paste0("./",c_resultsdirectory,"/model3/OUOU.model3.BestModel.RData"))
load(paste0("./",c_resultsdirectory,"/model4/OUOU.model4.BestModel.RData"))
load(paste0("./",c_resultsdirectory,"/model5/OUOU.model5.BestModel.RData"))
load(paste0("./",c_resultsdirectory,"/model6/OUOU.model6.BestModel.RData"))


best_model<-f_findBestModel(list(best_model1$best_output,best_model2$best_output,best_model3$best_output,best_model4$best_output,best_model5$best_output,best_model6$best_output))

save(best_model,file=paste0("./",c_resultsdirectory,"/BestFoundModel.RData"))

best_startpointforboot<-switch(
    best_model$best_index,
    best_startpointforboot_model1,
    best_startpointforboot_model2,
    best_startpointforboot_model3,
    best_startpointforboot_model4,
    best_startpointforboot_model5,
    best_startpointforboot_model6
)

best_parameter_signs<-switch(
    best_model$best_index,
    parameter_signs_model1,
    parameter_signs_model2,
    parameter_signs_model3,
    parameter_signs_model4,
    parameter_signs_model5,
    parameter_signs_model6
)    


best_Atype<-switch(
    best_model$best_index,
    Atype_model1,
    Atype_model2,
    Atype_model3,
    Atype_model4,
    Atype_model5,
    Atype_model6
)    

best_diagA<-switch(
    best_model$best_index,
    diagA_model1,
    diagA_model2,
    diagA_model3,
    diagA_model4,
    diagA_model5,
    diagA_model6
)   


mSyy<-best_model$best_output$ParamsInModel$Syy
best_SyyType<-"UpperTri"
if (sum(abs(mSyy[upper.tri(mSyy,diag=FALSE)])+abs(mSyy[lower.tri(mSyy,diag=FALSE)]))==0){best_SyyType<-"Diagonal"}

## Run 100 bootstrap repeats from best starting point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (run.parallel == TRUE) registerDoParallel(cores=num_cores)

OUOU.bootstrap.BestStart <- foreach(i=1:num_bootstrap_repeats) %dopar% {
    result <- multiResultClass()
    rexp(1) ## initialize .Random.seed  
    if (b_use_random_seed_from_manuscript){
        load(paste0("./",c_rseedoutfiledirectoryprefix,"/bootstrap/",c_rseedoutfiledirectoryprefix,"_bootstrap_startbest_iter_",i,".RData"))
        RNGkind(kind = RNG_kind_bootstrap_startbest_i[1], normal.kind = RNG_kind_bootstrap_startbest_i[2], sample.kind = RNG_kind_bootstrap_startbest_i[3]) 
        RNGversion(min(as.character(getRversion(),RNG_version_bootstrap_startbest_i)))
        rm(.Random.seed)
        assign('.Random.seed', random_seed_bootstrap_startbest_i, pos=.GlobalEnv)
    }else{
        RNG_kind_bootstrap_startbest_i<-RNGkind()
        RNG_version_bootstrap_startbest_i<-getRversion()
        random_seed_bootstrap_startbest_i<-.Random.seed
    }
    if (b_save_random_seeds){save(random_seed_bootstrap_startbest_i,RNG_kind_bootstrap_startbest_i,RNG_version_bootstrap_startbest_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/bootstrap/",c_rseedoutfiledirectoryprefix,"_bootstrap_startbest_iter_",i,".RData"))}
    start_time <- Sys.time()


  result$mvSLOUCH_bootstrap_output <- mvSLOUCH::parametric.bootstrap(estimated.model = best_model$best_output, phyltree = mvStree, 
                                                     regimes = regimes_Ellenberg,
                                                     Atype = best_Atype, 
                                                     Syytype = best_SyyType, 
                                                     diagA = best_diagA,
                                                     values.to.bootstrap = v_statistics_to_bootstrap,
                                                     numboot = 1,
                                                     start_point_for_optim=best_startpointforboot, 
                                                     parameter_signs=best_parameter_signs) 
    end_time <- Sys.time(); end_time - start_time 
    result$runtime <- end_time - start_time 
    res_iter_i<-result$mvSLOUCH_bootstrap_output;runtime<-result$runtime;save(res_iter_i,runtime,i,file=paste0("./",c_indivresultsdirectory ,"/bootstrap/OUOU.bootstrap.BestStart_iter_",i,".RData"))
    return(result)  
}
saveRDS(OUOU.bootstrap.BestStart, paste0("./",c_resultsdirectory,"/bootstrap/OUOU.bootstrap.BestStart.rds"))
## ==========  end of run 100 bootstrap repeats from best starting point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OUOU.bootstrap.RandomStart<-NA
if (!is.null(best_startpointforboot)){
    ## Run 100 bootstrap repeats from random starting point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (run.parallel == TRUE) registerDoParallel(cores=num_cores)

    OUOU.bootstrap.RandomStart <- foreach(i=1:num_bootstrap_repeats) %dopar% {
        result <- multiResultClass()
	rexp(1) ## initialize .Random.seed  
        if (b_use_random_seed_from_manuscript){
            load(paste0("./",c_rseedoutfiledirectoryprefix,"/bootstrap/",c_rseedoutfiledirectoryprefix,"_bootstrap_startrandom_iter_",i,".RData"))
	    RNGkind(kind = RNG_kind_bootstrap_startrandom_i[1], normal.kind = RNG_kind_bootstrap_startrandom_i[2], sample.kind = RNG_kind_bootstrap_startrandom_i[3]) 
            RNGversion(min(as.character(getRversion(),RNG_version_bootstrap_startrandom_i)))
            rm(.Random.seed)
            assign('.Random.seed', random_seed_bootstrap_startrandom_i, pos=.GlobalEnv)
        }else{
	    RNG_kind_bootstrap_startrandom_i<-RNGkind()
            RNG_version_bootstrap_startrandom_i<-getRversion()
	    random_seed_bootstrap_startrandom_i<-.Random.seed
        }
	if (b_save_random_seeds){save(random_seed_bootstrap_startrandom_i,RNG_kind_bootstrap_startrandom_i,RNG_version_bootstrap_startrandom_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/bootstrap/",c_rseedoutfiledirectoryprefix,"_bootstrap_startrandom_iter_",i,".RData"))}
        start_time <- Sys.time()


      result$mvSLOUCH_bootstrap_output <- mvSLOUCH::parametric.bootstrap(estimated.model = best_model$best_output, phyltree = mvStree, 
                                                     regimes = regimes_Ellenberg,
                                                     Atype = best_Atype, 
                                                     Syytype = best_SyyType, 
                                                     diagA = best_diagA,
                                                     values.to.bootstrap = v_statistics_to_bootstrap,
                                                     numboot = 1,
                                                     start_point_for_optim=NULL, 
                                                     parameter_signs=best_parameter_signs) 
	end_time <- Sys.time(); end_time - start_time 
        result$runtime <- end_time - start_time 
        res_iter_i<-result$mvSLOUCH_bootstrap_output;runtime<-result$runtime;save(res_iter_i,runtime,i,file=paste0("./",c_indivresultsdirectory ,"/bootstrap/OUOU.bootstrap.RandomStart_iter_",i,".RData"))
	return(result)  
    }
    saveRDS(OUOU.bootstrap.RandomStart, paste0("./",c_resultsdirectory,"/bootstrap/OUOU.bootstrap.RandomStart.rds"))
## ==========  end of run 100 bootstrap repeats from best starting point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


## ============= Extract bootstrapped statistics ==============================================================

f_get_bootci<-function(l_bootstrapped,name_stat){
    l_boottmp<-sapply(l_bootstrapped,function(x,name_stat){x$mvSLOUCH_bootstrap_output$bootstrapped.parameters[[which(names(x$mvSLOUCH_bootstrap_output$bootstrapped.parameters)==name_stat)]][[1]]},name_stat=name_stat,simplify=FALSE)
    m_tmp_LCI<-l_boottmp[[1]]
    if (is.matrix(m_tmp_LCI)){m_tmp_LCI[,]<-NA}
    else{
	if (is.list(m_tmp_LCI)){m_tmp_LCI<-sapply(m_tmp_LCI,function(x){x[]<-NA;x},simplify=FALSE)}
    }
    m_tmp_UCI<-m_tmp_LCI
    if (is.matrix(m_tmp_LCI)){
	for (i in 1:nrow(m_tmp_LCI)){
    	    for (j in 1:ncol(m_tmp_LCI)){
        	v_tmp<-sapply(l_boottmp,function(x,i,j){x[i,j]},i=i,j=j,simplify=TRUE)    
        	m_tmp_LCI[i,j]<-quantile(v_tmp,(1-bootci_lvl)/2)
        	m_tmp_UCI[i,j]<-quantile(v_tmp,bootci_lvl+(1-bootci_lvl)/2)
    	    }
	}
    }else{
        if (is.list(m_tmp_LCI)){
    	    for (i in 1:length(m_tmp_LCI)){
    		for (j in 1:length(m_tmp_LCI[[i]])){
        	    v_tmp<-sapply(l_boottmp,function(x,i,j){x[[i]][[j]]},i=i,j=j,simplify=TRUE)    
        	    m_tmp_LCI[[i]][[j]]<-quantile(v_tmp,(1-bootci_lvl)/2)
        	    m_tmp_UCI[[i]][[j]]<-quantile(v_tmp,bootci_lvl+(1-bootci_lvl)/2)
    		}
	    }
	}	
    }
    list(Lower_bootCI=m_tmp_LCI,param_est=NA,Upper_bootCI=m_tmp_UCI)
}


l_bootstrapped_stats_beststart<-vector("list",length(v_statistics_to_bootstrap))
l_bootstrapped_stats_randomstart<-vector("list",length(v_statistics_to_bootstrap))
names(l_bootstrapped_stats_beststart)<-v_statistics_to_bootstrap
names(l_bootstrapped_stats_randomstart)<-v_statistics_to_bootstrap


for (stat_to_boot in v_statistics_to_bootstrap){    
    l_bootstrapped_stats_beststart[[which(names(l_bootstrapped_stats_beststart)==stat_to_boot)]]<-f_get_bootci(OUOU.bootstrap.BestStart,stat_to_boot)
    l_bootstrapped_stats_beststart[[which(names(l_bootstrapped_stats_beststart)==stat_to_boot)]]$param_est<-best_model$best_output$ParamSummary[[which(names(best_model$best_output$ParamSummary)==stat_to_boot)]]
    if (!is.null(best_startpointforboot)){
	    l_bootstrapped_stats_randomstart[[which(names(l_bootstrapped_stats_randomstart)==stat_to_boot)]]<-f_get_bootci(OUOU.bootstrap.RandomStart,stat_to_boot)
	    l_bootstrapped_stats_randomstart[[which(names(l_bootstrapped_stats_randomstart)==stat_to_boot)]]$param_est<-best_model$best_output$ParamSummary[[which(names(best_model$best_output$ParamSummary)==stat_to_boot)]]
    }    
}
saveRDS(l_bootstrapped_stats_beststart, paste0("./",c_resultsdirectory,"/bootstrap/OUOU.bootstrapped_stats.BestStart.rds"))
if (!is.null(best_startpointforboot)){
    saveRDS(l_bootstrapped_stats_randomstart, paste0("./",c_resultsdirectory,"/bootstrap/OUOU.bootstrapped_stats.RandomStart.rds"))
}
## ============= End extract bootstrapped statistics ==========================================================



