## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

f_getres<-function(filename,simmodelist,b_returnbestmodels=FALSE){    
    load(filename)
    if (is.element("model_abbrev",names(simmodelist))){v_model_abbrev<-simmodelist$model_abbrev}
    if (is.element("model_setups",names(simmodelist))){simmodelist<-simmodelist$model_setups}    
    nummodsim<-length(simmodelist)
    if (b_returnbestmodels){
	lbestmodels<-vector("list",length(lres))
	lbestmodels<-sapply(lbestmodels,function(x,nummodsim){vector("list",nummodsim)},nummodsim=nummodsim,simplify=FALSE)
    }else{mResSummary<-matrix(0,nummodsim,nummodsim)}    
    k<-1
    for (res in lres){
	tmpmodeAICc<-rep(Inf,nummodsim)
	##if (b_returnbestmodels){ibestindex<-rep(NA,nummodsim)}
	j<-1
	for (tm in res$mvsl_estim_res$testedModels){
	    for (i in 1:length(simmodelist)){
		if (f_checksamemodel(simmodelist[[i]],tm$model)){
		    if (tm$aic.c<tmpmodeAICc[i]){
			tmpmodeAICc[i]<-tm$aic.c
			if (b_returnbestmodels){lbestmodels[[k]][[i]]<-tm;}
		    }
		}
	    }
	    j<-j+1
	}	
	if (!b_returnbestmodels && (length(simmodelist)>1)){	
	    for (i in 1:(length(simmodelist)-1)){
		for (j in (i+1):(length(simmodelist))){
		## it can happen that for a given estimation pair AICc will be the same for two models, these are excluded fromthe
		## counts as we only care about better-identifiability
		    if (tmpmodeAICc[i]<tmpmodeAICc[j]){
			mResSummary[i,j]<-mResSummary[i,j]+1
		    }
		    if (tmpmodeAICc[i]>tmpmodeAICc[j]){
			mResSummary[j,i]<-mResSummary[j,i]+1
		    }
		}
	    }
	}	
	k<-k+1
    }
    returnres<-NA
    if (b_returnbestmodels){
	returnres<-lbestmodels
    }else{
        rownames(mResSummary)<-v_model_abbrev #simmodelist$model_abbrev##recover_under$model_abbrev
        colnames(mResSummary)<-v_model_abbrev #simmodelist$model_abbrev##recover_under$model_abbrev
	returnres<-list(mResSummary=mResSummary,reps=length(lres))
    }
    returnres
}


f_checksamemodel<-function(model1,model2){
    b_isSame<-TRUE
    for (modfield in names(model1)){
	if (b_isSame){
	    if (is.element(modfield,names(model2))){
		i1<-which(names(model1)==modfield)
		i2<-which(names(model2)==modfield)
		if (is.null(model1[[i1]]) && (!is.null(model1[[i1]]))){b_isSame<-FALSE}
		else{
		    if ((!is.null(model1[[i1]])) && (is.null(model1[[i1]]))){b_isSame<-FALSE}
		    if (!((is.null(model1[[i1]])) && (is.null(model1[[i1]])))){
			if ((is.null(model2[[i2]]))||(model1[[i1]]!=model2[[i2]])){b_isSame<-FALSE}
		    }
		}
	    }else{b_isSame<-FALSE}
	}
    }
    b_isSame
}
