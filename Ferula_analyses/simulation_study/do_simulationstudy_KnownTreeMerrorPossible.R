## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

parSapply(cl,lmodelcomparison_setups,function(lsetup,vsetupstorun,runnum,main_directory,b_should_random_seed_be_saved,b_use_random_seed_from_manuscript){
#sapply(lmodelcomparison_setups,function(lsetup,vsetupstorun,runnum,main_directory,b_should_random_seed_be_saved,b_use_random_seed_from_manuscript){
if (vsetupstorun[[lsetup$id]]){
	source("tosource_for_mvSLOUCH.R")
	setuprun_directory<-paste0(main_directory,lsetup$id,"_ModelSelection/")
	dir.create(setuprun_directory,showWarnings=FALSE)
	
	sapply(lsetup$simulsetup$N,function(N,lsetup,runnum,main_directory,setuprun_directory){
	    ape_phyltree<-NA
	    lres<-NA
	    setup_comments_file<-paste0(setuprun_directory,"comments",runnum,"_",lsetup$id,".txt")
	    setup_error_comments_file<-paste0(setuprun_directory,"error_comments",runnum,"_",lsetup$id,".txt")
	    setuprun_directory_tmp<-paste0(setuprun_directory,"individual_runs_results_seeds/")
	    dir.create(setuprun_directory_tmp,showWarnings=FALSE)

	    v_iter<-lsetup$simulsetup$iter[[which(names(lsetup$simulsetup$iter)==paste0("size_",N))]]
	    
	    if (v_iter[1]>0){
		crandomseed_filename<-paste0(setuprun_directory,"RandomSeed_SetupID_",lsetup$id,"_N_",N,"_",runnum,".RData")
		rexp(1)
		if(b_use_random_seed_from_manuscript){##setup the random number generator if replicating results
		    b_randomseed_file_exists<-FALSE
		    if(file.exists(crandomseed_filename)){		    
			load(crandomseed_filename)
    			Rseed<-.Random.seed
			rm(.Random.seed)
    			RNGkind(kind = RNG_kind[1], normal.kind = RNG_kind[2], sample.kind = RNG_kind[3])
    			RNGversion(RNG_version)     
    			#.Random.seed<-Rseed
            		assign('.Random.seed', Rseed, pos=.GlobalEnv)
    			b_randomseed_file_exists<-TRUE
		    }
		}
		if((!b_use_random_seed_from_manuscript)||(!b_randomseed_file_exists)){
		    if (b_should_random_seed_be_saved){
			RNG_kind<-RNGkind()
			RNG_version<-getRversion()
			save(.Random.seed,RNG_kind,RNG_version,file=crandomseed_filename)
		    }
		}	
		if (v_iter[1]>1){## read in previous ones
		    lres_tmp<-vector("list",v_iter[1]-1)
		    for (j in 1:(v_iter[1]-1)){
			load(file=paste0(setuprun_directory_tmp,"Single_iteration_result_SetupID_",lsetup$id,"_N_",N,"_SimulationRun",runnum,"_",j,".RData"))
			lres_tmp[[j]]<-res_i
		    }
		}

		lres<-sapply(v_iter,function(i,N,lsetup,setup_comments_file,setup_error_comments_file,setuprun_directory,setuprun_directory_tmp,b_use_random_seed_from_manuscript,b_should_random_seed_be_saved){
        	    b_iter_not_done<-TRUE
		    attempt_counter<-1
		    while(b_iter_not_done){	
			clocalrandomseed_filename<-paste0(setuprun_directory_tmp,"RandomSeed_SetupID_",lsetup$id,"_N_",N,"_SimulationRun",runnum,"_",i,".RData")
			if(b_use_random_seed_from_manuscript){##setup the random number generator if replicating results
			    b_randomseed_file_exists<-FALSE

			    if(file.exists(clocalrandomseed_filename)){		    
				load(clocalrandomseed_filename)
				Rseed<-.Random.seed
    				rm(.Random.seed)
    				RNGkind(kind = RNG_kind[1], normal.kind = RNG_kind[2], sample.kind = RNG_kind[3])
    				RNGversion(RNG_version)     
    				#.Random.seed<-Rseed
    				assign('.Random.seed', Rseed, pos=.GlobalEnv)
    			        b_randomseed_file_exists<-TRUE
			    }			
			}
			if((!b_use_random_seed_from_manuscript)||(!b_randomseed_file_exists)){
			    if (b_should_random_seed_be_saved){
				RNG_kind<-RNGkind()
				RNG_version<-getRversion()
				save(.Random.seed,RNG_kind,RNG_version,file=clocalrandomseed_filename)
			    }
			}		
			tree_height<-1
			## monitoring of running time of process taken from https://stackoverflow.com/questions/34346619/how-to-stop-a-function-in-r-that-is-taking-too-long-and-give-it-an-alternative	
			base_time_limit<-5*60*60/1024
			time_limit <- 2*N*base_time_limit+600 ## time limit is in seconds
			setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)
			on.exit({setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)})
	    		    		    
			res_i<-NA
                        ape_phyltree<-NA
                        mData<-NA
                        mvsl_estim_res<-NA
                        lextracted_estim_summary<-NA
                        true_param_summary<-NA
                        Rseed_to_save_if_error<-NA
                        Merror<-NULL
                        tryCatch({
                    	    if (is.null(lsetup$ape_phyltree)){
                    		ape_phyltree<-TreeSim::sim.bd.taxa(N, 1, lambda=1, mu=0)[[1]]
                    		## rescale 
        			ape_phyltree$edge.length<-tree_height*ape_phyltree$edge.length/max(ape::node.depth.edgelength(ape_phyltree))
        		    }        	
        		    if (!is.null(lsetup$simulsetup$Merror)){
        			if (inherits(lsetup$simulsetup$Merror,"function")){
				    if (!is.null(lsetup$simulsetup$Merror_parameters)){
					Merror<-lsetup$simulsetup$Merror(N,lsetup$simulsetup$Merror_parameters)
				    }else{Merror<-lsetup$simulsetup$Merror(N)}
        			}else{Merror<-lsetup$simulsetup$Merror}
        		    }
    			    mData<-switch(lsetup$simulsetup$evolmodel,
				    bm=mvSLOUCH::simulBMProcPhylTree(ape_phyltree,X0=lsetup$simulsetup$modelparams$vX0,Sigma=lsetup$simulsetup$modelparams$Sxx,M.error=Merror),
    				    ouch=mvSLOUCH::simulOUCHProcPhylTree(ape_phyltree,lsetup$simulsetup$modelparams,M.error=Merror),
    				    mvslouch=mvSLOUCH::simulMVSLOUCHProcPhylTree(ape_phyltree,lsetup$simulsetup$modelparams,M.error=Merror)
    			    )
			    sink(setup_comments_file,append=TRUE);cat("Running simulation reestimation number: ");cat(i);cat(" for setup: ");cat(lsetup$id);cat(" and N=");cat(N);cat(" ");cat(format(Sys.time(), "%T %F"));cat("\n");sink()	
			    s_time<-proc.time()
                            Rseed_to_save_if_error<-.Random.seed

	    		    mvsl_estim_res<-mvSLOUCH::estimate.evolutionary.model(ape_phyltree, mData, regimes = NULL, 
    	    			    root.regime = NULL,  repeats =  lsetup$recoverunder$repeats, model.setups = lsetup$recoverunder$model_setups, 
        			    predictors = NULL, kY = lsetup$recoverunder$kY, doPrint = FALSE, pESS=NULL, 
            			    estimate.root.state=FALSE, min_bl = 0.0003, maxiter=c(10,50,100),M.error=Merror
            		    )
			    runtime<-proc.time()-s_time
			    lextracted_estim_summary<-extract_estimation_summary(lsetup$simulsetup,lsetup$recoverunder$model_setups,mvsl_estim_res)
	    		    true_param_summary<-switch(lsetup$simulsetup$evolmodel,
		    		    bm=mvSLOUCH::SummarizeBM(phyltree=ape_phyltree, mData=mData, modelParams=lextracted_estim_summary$TruePoint_model, t = c(1),M.error=Merror),
				    ouch=mvSLOUCH::SummarizeOUCH(phyltree=ape_phyltree, mData=mData, modelParams=lextracted_estim_summary$TruePoint_model, t = c(1),M.error=Merror),	
				    mvslouch=mvSLOUCH::SummarizeMVSLOUCH(phyltree=ape_phyltree, mData=mData, modelParams=lextracted_estim_summary$TruePoint_model, t = c(1),M.error=Merror)
			    )
			    b_iter_not_done<-FALSE
			}, error = function(e) {			
			    ## If there was some error with running the given iteration report it
			    sink(setup_error_comments_file,append=TRUE);cat("Error in simulation reestimation number: ");cat(i);cat(" for setup: ");cat(lsetup$id);cat(" and N=");cat(N);cat(" ");cat(format(Sys.time(), "%T %F "));cat(e$message);cat(" at attempt ");cat(attempt_counter);cat("\n");sink()	
			    b_iter_not_done<<-TRUE
		    	    attempt_counter<<-attempt_counter+1
		        })
		    }
		    res_i<-list(ape_phyltree=ape_phyltree,mData=mData,mvsl_estim_res=mvsl_estim_res,lextracted_estim_summary=lextracted_estim_summary,true_param_summary=true_param_summary,runtime=runtime,i=i,Merror=Merror)
		    save(res_i,file=paste0(setuprun_directory_tmp,"Single_iteration_result_SetupID_",lsetup$id,"_N_",N,"_SimulationRun",runnum,"_",i,".RData"))
		    res_i
		},N=N,lsetup=lsetup,setup_comments_file=setup_comments_file,setup_error_comments_file=setup_error_comments_file,setuprun_directory=setuprun_directory,setuprun_directory_tmp=setuprun_directory_tmp,b_use_random_seed_from_manuscript=b_use_random_seed_from_manuscript,b_should_random_seed_be_saved=b_should_random_seed_be_saved,simplify=FALSE)	

		if (v_iter[1]>1){## read in previous ones
    		    lres<-c(lres_tmp,lres)
		}
		mTrueChosen<-apply((t(sapply(lres,function(lrun){lrun$lextracted_estim_summary$TrueChosen},simplify=TRUE))),c(1,2),as.integer)
		vTrueChosenFrac<-apply(mTrueChosen,2,sum)/nrow(mTrueChosen)
	    
		save_file<-paste0(setuprun_directory,"SimulationReestimation_SetupID_",lsetup$id,"_N_",N,"_",runnum,".RData")
		save(lres,mTrueChosen,vTrueChosenFrac,file=save_file)
	    }
	    NA
	},lsetup=lsetup,runnum=runnum,main_directory=main_directory,setuprun_directory=setuprun_directory,simplify=FALSE)
    }
    NA   
},vsetupstorun=vsetupstorun,runnum=runnum,main_directory=main_directory,b_should_random_seed_be_saved=b_should_random_seed_be_saved,b_use_random_seed_from_manuscript=b_use_random_seed_from_manuscript,simplify=FALSE)



