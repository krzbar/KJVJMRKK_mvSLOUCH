## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


## MODEL 6  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Variable order is "plant.height", "seed.mass", "leaf.area", "leaf.mass"
parameter_signs_model6<-NULL
Atype_model6<-"DecomposablePositive"
diagA_model6<-NULL
## Run 100 Model 1 Diag's from the 100 OUOUstarts Diags. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (run.parallel == TRUE) registerDoParallel(cores=num_cores)
OUOUstart.Diags<-readRDS(paste0("./",c_resultsdirectory,"/initialdiagUTASR/OUOUstart.Diags.rds"))
OUOU.model6.Diags <- foreach(i=1:num_model_repeats) %dopar% {
    result <- multiResultClass()
    rexp(1) ## initialize .Random.seed  
    if (b_use_random_seed_from_manuscript){
        load(paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6diag_startdiag_iter_",i,".RData"))
        RNGkind(kind = RNG_kind_model6diag_startdiag_i[1], normal.kind = RNG_kind_model6diag_startdiag_i[2], sample.kind = RNG_kind_model6diag_startdiag_i[3]) 
        RNGversion(min(as.character(getRversion(),RNG_version_model6diag_startdiag_i)))
        rm(.Random.seed)
        assign('.Random.seed', random_seed_model6diag_startdiag_i, pos=.GlobalEnv)
    }else{
        RNG_kind_model6diag_startdiag_i<-RNGkind()
        RNG_version_model6diag_startdiag_i<-getRversion()
        random_seed_model6diag_startdiag_i<-.Random.seed
    }
    if (b_save_random_seeds){save(random_seed_model6diag_startdiag_i,RNG_kind_model6diag_startdiag_i,RNG_version_model6diag_startdiag_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6diag_startdiag_iter_",i,".RData"))}
    OUOUstart.Diags.i <- OUOUstart.Diags[[i]]$mvSLOUCH_output
    start_time <- Sys.time()
    result$mvSLOUCH_output <- mvSLOUCH::ouchModel(mvStree, mvData2, regimes=regimes_Ellenberg, root.regime=root_regime, 
                                                    Atype=Atype_model6, 
                                                    Syytype="Diagonal", 
                                                    diagA=diagA_model6, 
                                                    maxiter=c(10, 100), 
                                                    start_point_for_optim=list(
                                                      A=OUOUstart.Diags.i$FinalFound$ParamsInModel$A, 
                                                      Syy=OUOUstart.Diags.i$FinalFound$ParamsInModel$Syy), 
						    parameter_signs=parameter_signs_model6)
    end_time <- Sys.time(); end_time - start_time 
    result$runtime <- end_time - start_time 
    result$random_seed <-random_seed_model6diag_startdiag_i
    res_iter_i<-result$mvSLOUCH_output;runtime<-result$runtime;save(res_iter_i, runtime,i,file=paste0("./",c_indivresultsdirectory ,"/model6/OUOU.model6.Diags_iter_",i,".RData"))
    return(result)
}
saveRDS(OUOU.model6.Diags,  paste0("./",c_resultsdirectory,"/model6/OUOU.model6.Diags.rds"))
## ==========  end of run 100 Model 1 Diag's from the 100 OUOUstarts Diags. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Run 100 Model 1 Diag's from random starting points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (run.parallel == TRUE) registerDoParallel(cores=num_cores)

OUOU.model6.Diags.RandomStart <- foreach(i=1:num_model_repeats) %dopar% {
    result <- multiResultClass()
    rexp(1) ## initialize .Random.seed  
    if (b_use_random_seed_from_manuscript){
        load(paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6diag_startrandom_iter_",i,".RData"))
        RNGkind(kind = RNG_kind_model6diag_startrandom_i[1], normal.kind = RNG_kind_model6diag_startrandom_i[2], sample.kind = RNG_kind_model6diag_startrandom_i[3]) 
        RNGversion(min(as.character(getRversion(),RNG_version_model6diag_startrandom_i)))
        rm(.Random.seed)
        assign('.Random.seed', random_seed_model6diag_startrandom_i, pos=.GlobalEnv)
    }else{
        RNG_kind_model6diag_startrandom_i<-RNGkind()
        RNG_version_model6diag_startrandom_i<-getRversion()
        random_seed_model6diag_startrandom_i<-.Random.seed
    }
    if (b_save_random_seeds){save(random_seed_model6diag_startrandom_i,RNG_kind_model6diag_startrandom_i,RNG_version_model6diag_startrandom_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6diag_startrandom_iter_",i,".RData"))}
    start_time <- Sys.time()
    result$mvSLOUCH_output <- mvSLOUCH::ouchModel(mvStree, mvData2, regimes=regimes_Ellenberg, root.regime=root_regime, 
                                                    Atype=Atype_model6, 
                                                    Syytype="Diagonal", 
                                                    diagA=diagA_model6, 
                                                    maxiter=c(10, 100), 
                                                    start_point_for_optim=NULL, 
                                                    parameter_signs=parameter_signs_model6)  
    end_time <- Sys.time(); end_time - start_time 
    result$runtime <- end_time - start_time 
    result$random_seed <-random_seed_model6diag_startrandom_i
    res_iter_i<-result$mvSLOUCH_output;runtime<-result$runtime;save(res_iter_i,runtime,i,file=paste0("./",c_indivresultsdirectory ,"/model6/OUOU.model6.Diags.RandomStart_iter_",i,".RData"))
    return(result)  
}
saveRDS(OUOU.model6.Diags.RandomStart,  paste0("./",c_resultsdirectory,"/model6/OUOU.model6.Diags.RandomStart.rds"))
## ==========  end of run 100 Model 1 Diag's from random starting points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ========= Extract best diagonal model for starting point of upper triangualr analyses =================================
best_diag_extracted<-f_extractBestFromPair(OUOU.model6.Diags,OUOU.model6.Diags.RandomStart)
best_diag_startingpoint<-best_diag_extracted$best_model_startingpoint
best_startpointforboot_diag<-NULL
if (best_diag_extracted$which_of_pair==1){best_startpointforboot_diag<-f_createstartingpoint(OUOUstart.Diags[[best_diag_extracted$best_index]]$mvSLOUCH_output)}
## ========= End of extract best diagonal model for starting point of upper triangualr analyses ==========================

## MODEL 6  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run 100 Model 6 UT's from the 100 OUOUstarts UT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (run.parallel == TRUE) registerDoParallel(cores=num_cores)
OUOUstart.UT<-readRDS(paste0("./",c_resultsdirectory,"/initialdiagUTASR/OUOUstart.UT.rds"))

OUOU.model6.UT <- foreach(i=1:num_model_repeats) %dopar% {
  
    result <- multiResultClass()
    rexp(1) ## initialize .Random.seed  
    if (b_use_random_seed_from_manuscript){
        load(paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6UT_startUT_iter_",i,".RData"))
        RNGkind(kind = RNG_kind_model6UT_startUT_i[1], normal.kind = RNG_kind_model6UT_startUT_i[2], sample.kind = RNG_kind_model6UT_startUT_i[3]) 
        RNGversion(min(as.character(getRversion(),RNG_version_model6UT_startUT_i)))
        rm(.Random.seed)
        assign('.Random.seed', random_seed_model6UT_startUT_i, pos=.GlobalEnv)
    }else{
        RNG_kind_model6UT_startUT_i<-RNGkind()
        RNG_version_model6UT_startUT_i<-getRversion()
        random_seed_model6UT_startUT_i<-.Random.seed
    }
    if (b_save_random_seeds){save(random_seed_model6UT_startUT_i,RNG_kind_model6UT_startUT_i,RNG_version_model6UT_startUT_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6UT_startUT_iter_",i,".RData"))}


    OUOUstart.UT.i <- OUOUstart.UT[[i]]$mvSLOUCH_output
    start_time <- Sys.time()
    result$mvSLOUCH_output <- mvSLOUCH::ouchModel(mvStree, mvData2, regimes=regimes_Ellenberg,
                                                Atype=Atype_model6, 
                                                Syytype="UpperTri", diagA=diagA_model6, 
                                                maxiter=c(10, 100), 
                                                start_point_for_optim=list(
                                                  A=OUOUstart.UT.i$FinalFound$ParamsInModel$A, 
                                                  Syy=OUOUstart.UT.i$FinalFound$ParamsInModel$Syy), 
                                                parameter_signs=parameter_signs_model6 )

    end_time <- Sys.time(); end_time - start_time 
    result$runtime <- end_time - start_time 
    result$random_seed <-random_seed_model6UT_startUT_i
    res_iter_i<-result$mvSLOUCH_output;runtime<-result$runtime;save(res_iter_i, runtime,i,file=paste0("./",c_indivresultsdirectory ,"/model6/OUOU.model6.UT_iter_",i,".RData"))
    return(result)  
}
saveRDS(OUOU.model6.UT, paste0("./",c_resultsdirectory,"/model6/OUOU.model6.UT.rds"))
## ==========  end of run 100 Model 6 UT's from the 100 OUOUstarts UT. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Run 100 Model 1 UT's from random starting points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (run.parallel == TRUE) registerDoParallel(cores=num_cores)

OUOU.model6.UT.RandomStart <- foreach(i=1:num_model_repeats) %dopar% {
    result <- multiResultClass()
    rexp(1) ## initialize .Random.seed  
    if (b_use_random_seed_from_manuscript){
        load(paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6UT_startrandom_iter_",i,".RData"))
        RNGkind(kind = RNG_kind_model6UT_startrandom_i[1], normal.kind = RNG_kind_model6UT_startrandom_i[2], sample.kind = RNG_kind_model6UT_startrandom_i[3]) 
        RNGversion(min(as.character(getRversion(),RNG_version_model6UT_startrandom_i)))
        rm(.Random.seed)
        assign('.Random.seed', random_seed_model6UT_startrandom_i, pos=.GlobalEnv)
    }else{
        RNG_kind_model6UT_startrandom_i<-RNGkind()
        RNG_version_model6UT_startrandom_i<-getRversion()
        random_seed_model6UT_startrandom_i<-.Random.seed
    }
    if (b_save_random_seeds){save(random_seed_model6UT_startrandom_i,RNG_kind_model6UT_startrandom_i,RNG_version_model6UT_startrandom_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6UT_startrandom_iter_",i,".RData"))}
    start_time <- Sys.time()

    result$mvSLOUCH_output <- mvSLOUCH::ouchModel(mvStree, mvData2, regimes=regimes_Ellenberg, root.regime=root_regime, 
                                                Atype=Atype_model6, 
                                                Syytype="UpperTri", diagA=diagA_model6, 
                                                maxiter=c(10, 100), 
                                                start_point_for_optim=NULL, 
                                                parameter_signs=parameter_signs_model6)
    end_time <- Sys.time(); end_time - start_time 
    result$runtime <- end_time - start_time 
    result$random_seed <-random_seed_model6UT_startrandom_i
    res_iter_i<-result$mvSLOUCH_output;runtime<-result$runtime;save(res_iter_i,runtime,i,file=paste0("./",c_indivresultsdirectory ,"/model6/OUOU.model6.UT.RandomStart_iter_",i,".RData"))
    return(result)  
}
saveRDS(OUOU.model6.UT.RandomStart, paste0("./",c_resultsdirectory,"/model6/OUOU.model6.UT.RandomStart.rds"))
## ==========  end of run 100 Model 6 UT from random starting points ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ========= Extract best UT model for starting point of upper triangualr analyses =================================
best_UT_extracted<-f_extractBestFromPair(OUOU.model6.UT,OUOU.model6.UT.RandomStart)
## ========= End of extract best diagonal model for starting point of upper triangualr analyses ==========================



## Run Model 6 UT's from best found diag starting point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OUOU.model6.UT.BestDiagStart <- list() 
    result <- multiResultClass()
    rexp(1) ## initialize .Random.seed  
    i<-0;
    if (b_use_random_seed_from_manuscript){
        load(paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6UT_startbestdiag_iter_",i,".RData"))
        RNGkind(kind = RNG_kind_model6UT_startbestdiag_i[1], normal.kind = RNG_kind_model6UT_startbestdiag_i[2], sample.kind = RNG_kind_model6UT_startbestdiag_i[3]) 
        RNGversion(min(as.character(getRversion(),RNG_version_model6UT_startbestdiag_i)))
        rm(.Random.seed)
        assign('.Random.seed', random_seed_model6UT_startbestdiag_i, pos=.GlobalEnv)
    }else{
        RNG_kind_model6UT_startbestdiag_i<-RNGkind()
        RNG_version_model6UT_startbestdiag_i<-getRversion()
        random_seed_model6UT_startbestdiag_i<-.Random.seed
    }
    if (b_save_random_seeds){save(random_seed_model6UT_startbestdiag_i,RNG_kind_model6UT_startbestdiag_i,RNG_version_model6UT_startbestdiag_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/model6/",c_rseedoutfiledirectoryprefix,"_model6UT_startbestdiag_iter_",i,".RData"))}
    start_time <- Sys.time()

    result$mvSLOUCH_output <- mvSLOUCH::ouchModel(mvStree, mvData2, regimes=regimes_Ellenberg, root.regime=root_regime, 
                                                Atype=Atype_model6, 
                                                Syytype="UpperTri", diagA=diagA_model6, 
                                                maxiter=c(10, 100), 
                                                start_point_for_optim=best_diag_startingpoint, 
                                                parameter_signs=parameter_signs_model6)
    end_time <- Sys.time(); end_time - start_time 
    result$runtime <- end_time - start_time 
    result$random_seed <-random_seed_model6UT_startbestdiag_i
    res_iter_i<-result$mvSLOUCH_output;runtime<-result$runtime;save(res_iter_i,runtime,i,file=paste0("./",c_indivresultsdirectory ,"/model6/OUOU.model6.UT.BestDiagStart_iter_",i,".RData"))

OUOU.model6.UT.BestDiagStart <- result
saveRDS(OUOU.model6.UT.BestDiagStart, paste0("./",c_resultsdirectory,"/model6/OUOU.model6.UT.BestDiagStart.rds"))
## ==========  end of run of Model 6 UT from best found diag starting point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ========= Extract best UT model ==========================================================================================
best_UT_extracted<-f_extractBestFromPair(OUOU.model6.UT,OUOU.model6.UT.RandomStart)
best_startpointforboot_UT<-NULL
if (best_UT_extracted$which_of_pair==1){best_startpointforboot_UT<-f_createstartingpoint(OUOUstart.UT[[best_UT_extracted$best_index]]$mvSLOUCH_output)}
## ==========================================================================================================================

## ========= Extract best model =============================================================================================
best_model6<-f_findBestModel(list(best_diag_extracted$best_model,best_UT_extracted$best_model,OUOU.model6.UT.BestDiagStart$mvSLOUCH_output))
best_startpointforboot_model6<-switch(
    best_model6$best_index,
    best_startpointforboot_diag,
    best_startpointforboot_UT,
    best_diag_startingpoint
)
## ==========================================================================================================================

save(best_model6,best_startpointforboot_model6,parameter_signs_model6,Atype_model6,diagA_model6,file=paste0("./",c_resultsdirectory,"/model6/OUOU.model6.BestModel.RData"))
    