## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## ====================== Run first initial model ========================================================
## 100 OUOUstart.Diag models to provide 100 starting points for the six Syytype="Diagonal" models

if (run.parallel == TRUE){registerDoParallel(cores=num_cores)}

if (b_use_random_seed_from_manuscript){.Random.seed<-random_seed_initial}else{random_seed_initial<-.Random.seed}
if (b_save_random_seeds){save(random_seed_initial,random_seed_ASR,RNG_kind,RNG_version,file=paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASR_initial.RData"))}

starting_diagA<-NULL
OUOUstart.Diags <- foreach(i=1:num_model_repeats) %dopar% {
    rexp(1) ## initialize .Random.seed
    result <- multiResultClass()
  
    if (b_use_random_seed_from_manuscript) {
	load(paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_initialdiag_iter_",i,".RData"))
        RNGkind(kind = RNG_kind_initialdiag_i[1], normal.kind = RNG_kind_initialdiag_i[2], sample.kind = RNG_kind_initialdiag_i[3]) 
	RNGversion(min(as.character(getRversion(),RNG_version_initialdiag_i)))
	rm(.Random.seed)
	assign('.Random.seed', random_seed_initialdiag_i, pos=.GlobalEnv)
    }else{
	RNG_kind_initialdiag_i<-RNGkind()
	RNG_version_initialdiag_i<-getRversion()
	random_seed_initialdiag_i<-.Random.seed
    }
    if (b_save_random_seeds){save(random_seed_initialdiag_i,RNG_kind_initialdiag_i,RNG_version_initialdiag_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_initialdiag_iter_",i,".RData"))}

    start_time <- Sys.time()
    result$mvSLOUCH_output <- mvSLOUCH::ouchModel( mvStree, mvData2, regimes=regimes_Ellenberg, root.regime=root_regime, 
                                         Atype="DecomposablePositive",
                                         Syytype="Diagonal", diagA=starting_diagA,
                                         maxiter=c(10, 100) )
    end_time <- Sys.time(); end_time - start_time 
    result$runtime <- end_time - start_time  
    result$random_seed <-random_seed_initialdiag_i
    res_iter_i<-result$mvSLOUCH_output;runtime<-result$runtime;save(res_iter_i, runtime,i,file=paste0("./",c_indivresultsdirectory ,"/initialdiagUTASR/OUOUstart.Diags_iter_",i,".RData"))
    return(result)  
}
saveRDS(OUOUstart.Diags, paste0("./",c_resultsdirectory,"/initialdiagUTASR/OUOUstart.Diags.rds"))

## 100 OUOUstart.UT models to provide 100 starting points for the six Syytype="UpperTri" models

if (run.parallel == TRUE){registerDoParallel(cores=num_cores)}

OUOUstart.UT <- foreach(i=1:num_model_repeats) %dopar% {
    rexp(1) ## initialize .Random.seed
    result <- multiResultClass()
  
    if (b_use_random_seed_from_manuscript) {
	load(paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_initialUT_iter_",i,".RData"))
        RNGkind(kind = RNG_kind_initialUT_i[1], normal.kind = RNG_kind_initialUT_i[2], sample.kind = RNG_kind_initialUT_i[3]) 
	RNGversion(min(as.character(getRversion(),RNG_version_initialUT_i)))	
	rm(.Random.seed)
	assign('.Random.seed', random_seed_initialUT_i, pos=.GlobalEnv)
    }else{
	RNG_kind_initialUT_i<-RNGkind()
	RNG_version_initialUT_i<-getRversion()
	random_seed_initialUT_i<-.Random.seed
    }
    if (b_save_random_seeds){save(random_seed_initialUT_i,RNG_kind_initialUT_i,RNG_version_initialUT_i,file=paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_initialUT_iter_",i,".RData"))}

    start_time <- Sys.time()

    result$mvSLOUCH_output <- mvSLOUCH::ouchModel( mvStree, mvData2, regimes=regimes_Ellenberg, root.regime=root_regime, 
                                         Atype="DecomposablePositive",
                                         Syytype="UpperTri", diagA=starting_diagA,
                                         maxiter=c(10, 100) )

    end_time <- Sys.time(); end_time - start_time 
    result$runtime <- end_time - start_time  
    result$random_seed <-random_seed_initialUT_i
    res_iter_i<-result$mvSLOUCH_output;runtime<-result$runtime;save(res_iter_i, runtime,i,file=paste0("./",c_indivresultsdirectory ,"/initialdiagUTASR/OUOUstart.UT_iter_",i,".RData"))
    return(result)  
}
saveRDS(OUOUstart.UT, paste0("./",c_resultsdirectory,"/initialdiagUTASR/OUOUstart.UT.rds"))
## ====================== End run first initial model =====================================================
