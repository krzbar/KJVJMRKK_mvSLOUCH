## This file accompanies the manuscript:  Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje "Analytical advances alleviate model misspecification in non--Brownian multivariate comparative methods"


## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## This script contains the functions needed to estimate under independent variables Brownian motion.
## Its functionality will be ported to mvSLOUCH in the future.

BrownianMotionModel_independent<-function(phyltree,mData,predictors=NULL,M.error=NULL,min_bl=0.0003){
    p<-ncol(mData)
    numobs<-p*nrow(mData)    
    numobs<- numobs-length(which(is.na(c(mData))))
    
    l_BM1d_res<-vector("list",p)
    v_trait_names<-NA
    if (!is.null(colnames(mData))){
	names(l_BM1d_res)<-colnames(mData)
	v_trait_names<-colnames(mData)
    }
    for (i in 1:p){
	mData_1D<-mData[,i,drop=FALSE]
	l_BM1d_res[[i]]<-BrownianMotionModel(phyltree,mData_1D,predictors=predictors,M.error=M.error,min_bl=min_bl)
    }

    BM_res<-list(ParamsInModel=list(vX0=matrix(0,ncol=1,nrow=p),Sxx=matrix(0,p,p)),ParamSummary=list("StS"=matrix(0,p,p),"LogLik"=0,"dof"=0,"m2loglik"=0,"aic"=0,"aic.c"=0,"sic"=0,"bic"=0,"RSS"=list(RSS=0,R2=0),"confidence.interval"=list(regression.summary=list(X0.regression.confidence.interval=matrix(0,ncol=3,nrow=p),regression.covariance.matrix=matrix(0,p,p),regression.confidence.interval.comment=""))))
    if (!is.na(v_trait_names[1])){
	rownames(BM_res$ParamsInModel$vX0)<-v_trait_names
	rownames(BM_res$ParamsInModel$Sxx)<-v_trait_names
	colnames(BM_res$ParamsInModel$Sxx)<-v_trait_names
	rownames(BM_res$ParamSummary$StS)<-v_trait_names
	colnames(BM_res$ParamSummary$StS)<-v_trait_names
    }

    colnames(BM_res$ParamSummary$confidence.interval$regression.summary$X0.regression.confidence.interval)<-c("Lower.end", "Estimated.Point", "Upper.end")
    rownames(BM_res$ParamSummary$confidence.interval$regression.summary$regression.covariance.matrix)<-paste0("X0_",1:p)
    colnames(BM_res$ParamSummary$confidence.interval$regression.summary$regression.covariance.matrix)<-paste0("X0_",1:p)
    vRSS0<-rep(Inf,p)
    for (i in 1:p){
	## the results
	BM_res$ParamsInModel$vX0[i]<-l_BM1d_res[[i]]$ParamsInModel$vX0[1,1]
	BM_res$ParamsInModel$Sxx[i,i]<-l_BM1d_res[[i]]$ParamsInModel$Sxx[1,1]
	
	## the summary, can be also extracted from the output of SummarizeBM
#	BM_res$ParamSummary$StS[i,i]<-l_BM1d_res[[i]]$ParamSummary$StS[1,1]	
#	BM_res$ParamSummary$LogLik<-BM_res$ParamSummary$LogLik+l_BM1d_res[[i]]$ParamSummary$LogLik
	BM_res$ParamSummary$dof<-BM_res$ParamSummary$dof+l_BM1d_res[[i]]$ParamSummary$dof
#	BM_res$ParamSummary$m2loglik<-BM_res$ParamSummary$m2loglik+l_BM1d_res[[i]]$ParamSummary$m2loglik
#	BM_res$ParamSummary$aic<-BM_res$ParamSummary$aic+l_BM1d_res[[i]]$ParamSummary$aic
	BM_res$ParamSummary$RSS$RSS<-BM_res$ParamSummary$RSS$RSS+l_BM1d_res[[i]]$ParamSummary$RSS$RSS
	vRSS0[i]<-l_BM1d_res[[i]]$ParamSummary$RSS$RSS/(1-l_BM1d_res[[i]]$ParamSummary$RSS$R2)
	BM_res$ParamSummary$confidence.interval$regression.summary$X0.regression.confidence.interval[i,]<-l_BM1d_res[[i]]$ParamSummary$confidence.interval$regression.summary$X0.regression.confidence.interval[1,]
	BM_res$ParamSummary$confidence.interval$regression.summary$regression.covariance.matrix[i,i]<-l_BM1d_res[[i]]$ParamSummary$confidence.interval$regression.summary$regression.covariance.matrix[1,1]
	BM_res$ParamSummary$confidence.interval$regression.summary$regression.confidence.interval.comment<-l_BM1d_res[[i]]$ParamSummary$confidence.interval$regression.summary$regression.confidence.interval.comment
    }
    BM_res$ParamSummary$RSS$R2<-1-BM_res$ParamSummary$RSS$RSS/(sum(vRSS0))
#    BM_res$ParamSummary$aic.c<- BM_res$ParamSummary$aic +2*BM_res$ParamSummary$dof*(BM_res$ParamSummary$dof+1)/(numobs-BM_res$ParamSummary$dof-1)
#    BM_res$ParamSummary$sic<- BM_res$ParamSummary$m2loglik+log(numobs)*BM_res$ParamSummary$dof
#    BM_res$ParamSummary$bic<-BM_res$ParamSummary$m2loglik+BM_res$ParamSummary$dof*log(numobs)
    ## extract the remaining data from SummarizeBM
    tmp_BMres<-mvSLOUCH::SummarizeBM(phyltree=phyltree, mData=mData, modelParams=BM_res$ParamsInModel,dof=BM_res$ParamSummary$dof)[[1]]
    BM_res$ParamSummary$StS<-tmp_BMres$PointSummary$StS
    BM_res$ParamSummary$LogLik<-tmp_BMres$PointSummary$LogLik
    BM_res$ParamSummary$m2loglik<-tmp_BMres$PointSummary$m2loglik
    BM_res$ParamSummary$aic<-tmp_BMres$PointSummary$aic
    BM_res$ParamSummary$aic.c<-tmp_BMres$PointSummary$aic.c
    BM_res$ParamSummary$sic<-tmp_BMres$PointSummary$sic
    BM_res$ParamSummary$bic<-tmp_BMres$PointSummary$bic
    BM_res
}

BrownianMotionModel_independent_estim_res<-function(phyltree,mData,M.error=NULL){
    BM_res<-BrownianMotionModel_independent(phyltree,mData,M.error)
    list(result=BM_res,aic.c=BM_res$ParamSummary$aic.c,bic=BM_res$ParamSummary$bic,model=list(evolmodel="bm_indep"))
}

join_estimation_results<-function(l_mvsl_estim_res){
## function has to be a list of outputs of mvSLOUCH::estimate.evolutionary.model()
    joined_mvsl_estim_res<-list(BestModel=list("model.description"=NULL,"key.properties"=NULL,"parameter.SE"=NULL,"aic.c"=Inf,"evolmodel"=NULL,"model"=NULL,"model.call"=NA,"BestModel"=NULL,"i"=0),testedModels=list(),model.setups=list(),repeats=0)
    num_res_lists<-length(l_mvsl_estim_res)
    joined_mvsl_estim_res$testedModels<-l_mvsl_estim_res[[1]]$testedModels
    joined_mvsl_estim_res$model.setups<-l_mvsl_estim_res[[1]]$model.setups
    joined_mvsl_estim_res$repeats<-l_mvsl_estim_res[[1]]$repeats


    if (num_res_lists>1){
        for (i in 2:num_res_lists){
    	    joined_mvsl_estim_res$testedModels<-c(joined_mvsl_estim_res$testedModels,l_mvsl_estim_res[[i]]$testedModels)
    	    joined_mvsl_estim_res$model.setups<-c(joined_mvsl_estim_res$model.setups,l_mvsl_estim_res[[i]]$model.setups)
    	    if (joined_mvsl_estim_res$repeats!=l_mvsl_estim_res[[i]]$repeats){
    		print("WARNING: number of repeats different between different analyses, taking minimum!")
    		joined_mvsl_estim_res$repeats<-min(joined_mvsl_estim_res$repeats,l_mvsl_estim_res[[i]]$repeats)
    	    }
	} 
    }
    base_i<-0
    for (i in 1:num_res_lists){
	if (is.list(l_mvsl_estim_res[[i]]$BestModel)&&(is.na(l_mvsl_estim_res[[i]]$BestModel)[1])){
	    if(l_mvsl_estim_res[[i]]$BestModel$aic.c<joined_mvsl_estim_res$BestModel$aic.c){
		joined_mvsl_estim_res$BestModel$aic.c<-l_mvsl_estim_res[[i]]$BestModel$aic.c
		joined_mvsl_estim_res$BestModel$i<-base_i+l_mvsl_estim_res[[i]]$BestModel$i
	    }
	}else{
	    for (j in 1:length(l_mvsl_estim_res[[i]]$testedModels)){
		if (is.list(l_mvsl_estim_res[[i]]$testedModels[[j]])&&(!is.na(l_mvsl_estim_res[[i]]$testedModels[[j]])[1])&&(!is.na(l_mvsl_estim_res[[i]]$testedModels[[j]]$aic.c))&&(!is.null(l_mvsl_estim_res[[i]]$testedModels[[j]]$aic.c))){
		    if(l_mvsl_estim_res[[i]]$testedModels[[j]]$aic.c<joined_mvsl_estim_res$BestModel$aic.c){
			joined_mvsl_estim_res$BestModel$aic.c<-l_mvsl_estim_res[[i]]$testedModels[[j]]$aic.c
			joined_mvsl_estim_res$BestModel$evolmodel<-l_mvsl_estim_res[[i]]$testedModels[[j]]$model$evolmodel
			joined_mvsl_estim_res$BestModel$BestModel<-l_mvsl_estim_res[[i]]$testedModels[[j]]$result
			joined_mvsl_estim_res$BestModel$model<-l_mvsl_estim_res[[i]]$testedModels[[j]]$model
			joined_mvsl_estim_res$BestModel$i<-base_i+j			
		    }
		}
	    }
	}	
        base_i<-base_i+length(l_mvsl_estim_res[[i]]$testedModels)
    }    
    joined_mvsl_estim_res
}

