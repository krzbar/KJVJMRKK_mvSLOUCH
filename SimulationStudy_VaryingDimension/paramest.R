## This file accompanies the manuscript:  Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje "Rotation Invariance in non--Brownian motion Phylogenetic Comparative Methods: A comment on Adams and Collyer (2018)"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

f_findvalueinlist<-function(lmodelparams,vlocation){
    ltruelist<-lmodelparams
    res<-NULL
    for (cfieldname in vlocation){
        if (cfieldname=="FinalFound"){
    	    if (is.list(ltruelist$MaxLikFound)){
    		cfieldname<-"MaxLikFound"
    	    }
        }
        v_which<-which(names(ltruelist)==cfieldname)
        if (length(v_which)>0){ltruelist<-ltruelist[[v_which]];res<-ltruelist}
        else{res<-NULL}        
    }
    res
}

## obtaining model parameters from setup file
vmodelids<-sapply(lmodelcomparison_setups,function(x){x$id},simplify=TRUE)
lsimmodellist<-sapply(lsimmodellist,function(x,lmodelcomparison_setups,vmodelids){
    model_index<-which(vmodelids==x$id)
    x$modelparams<-lmodelcomparison_setups[[model_index]]$simulsetup$modelparams
    x
},lmodelcomparison_setups=lmodelcomparison_setups,vmodelids=vmodelids,simplify=FALSE)



## calculate the relative distance between two vectors, matrices
## first argument, x is the true value
## if it is 0, then just the Euclidean distance is returned
fdistL2scalarproduct<-function(x,y){sqrt(((x-y)%*%(x-y))[1,1])}
fdistL2scalarproduct_relative<-function(x,y){x2<-(x%*%x)[1,1];if (isTRUE(all.equal(x2,0))){res<-fdistL2scalarproduct(x,y)}else{res<-sqrt(((x-y)%*%(x-y))[1,1])/sqrt(x2)};res}
fL2<-function(x,y){x<-c(x);y<-c(y);fdistL2scalarproduct(x,y)}
fL2list<-function(x,y){x<-unlist(x);y<-unlist(y);fdistL2scalarproduct(x,y)}
fL2_relative<-function(x,y){#print(x);print(y);print("====");
x<-c(x);y<-c(y);fdistL2scalarproduct_relative(x,y)}
fL2list_relative<-function(x,y){x<-unlist(x);y<-unlist(y);fdistL2scalarproduct_relative(x,y)}
fL2_used<-fL2_relative
fL2list_used<-fL2list_relative



lparamstocheck<-list(
list(param="vX0",evolmodel="bm",location_true=c("modelparams","vX0"),location_est=c("ParamsInModel","vX0"),nameforplot="X0",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,3),xlab="",ylab="",main="X(0)")
,list(param="Sxx",evolmodel="bm",location_true=c("modelparams","Sxx"),location_est=c("ParamsInModel","Sxx"),nameforplot="Sxx",f_cf=fL2_used,paramtype="",ylim=c(0,10),xlim=c(0,0.2),xlab="",ylab="",main=expression(Sigma[xx]))
,list(param="vX0",evolmodel="bm_indep",location_true=c("modelparams","vX0"),location_est=c("ParamsInModel","vX0"),nameforplot="X0",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,3),xlab="",ylab="",main="X(0)")
,list(param="Sxx",evolmodel="bm_indep",location_true=c("modelparams","Sxx"),location_est=c("ParamsInModel","Sxx"),nameforplot="Sxx",f_cf=fL2_used,paramtype="",ylim=c(0,10),xlim=c(0,0.2),xlab="",ylab="",main=expression(Sigma[xx]))
,list(param="A",evolmodel="ouch",location_true=c("modelparams","A"),location_est=c("FinalFound","ParamsInModel","A"),nameforplot="A",f_cf=fL2_used,paramtype="",ylim=c(0,1),xlim=c(0,2),xlab="",ylab="",main="A")
,list(param="Syy",evolmodel="ouch",location_true=c("modelparams","Syy"),location_est=c("FinalFound","ParamsInModel","Syy"),nameforplot="Syy",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,0.1),xlab="",ylab="",main=expression(Sigma[yy]))
,list(param="vY0",evolmodel="ouch",location_true=c("modelparams","vY0"),location_est=c("FinalFound","ParamsInModel","vY0"),nameforplot="Y0",f_cf=fL2_used,paramtype="",ylim=c(0,1),xlim=c(0,10),xlab="",ylab="",main="Y(0)")
,list(param="vY0",evolmodel="mvslouch",location_true=c("modelparams","vY0"),location_est=c("FinalFound","ParamsInModel","vY0"),nameforplot="Y0_BML",f_cf=fL2_used,paramtype="",ylim=c(0,1),xlim=c(0,10),xlab="",ylab="",main="Y(0)")
,list(param="mPsi",evolmodel="ouch",location_true=c("modelparams","mPsi"),location_est=c("FinalFound","ParamsInModel","mPsi"),nameforplot="Psi",f_cf=fL2_used,paramtype="",ylim=c(0,0.5),xlim=c(0,10),xlab="",ylab="",main=expression(Psi))
,list(param="mPsi",evolmodel="mvslouch",location_true=c("modelparams","mPsi"),location_est=c("FinalFound","ParamsInModel","mPsi"),nameforplot="Psi_BML",f_cf=fL2_used,paramtype="",ylim=c(0,0.5),xlim=c(0,10),xlab="",ylab="",main=expression(Psi))
,list(param="A",evolmodel="mvslouch",location_true=c("modelparams","A"),location_est=c("FinalFound","ParamsInModel","A"),nameforplot="A_BML",f_cf=fL2_used,paramtype="",ylim=c(0,0.5),xlim=c(0,20),xlab="",ylab="",main="A")
,list(param="Syy",evolmodel="mvslouch",location_true=c("modelparams","Syy"),location_est=c("FinalFound","ParamsInModel","Syy"),nameforplot="Syy_BML",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,0.1),xlab="",ylab="",main=expression(Sigma[yy]))
,list(param="Sxx",evolmodel="mvslouch",location_true=c("modelparams","Sxx"),location_est=c("FinalFound","ParamsInModel","Sxx"),nameforplot="Sxx_BML",f_cf=fL2_used,paramtype="",ylim=c(0,10),xlim=c(0,0.2),xlab="",ylab="",main=expression(Sigma[xx]))
,list(param="vX0",evolmodel="mvslouch",location_true=c("modelparams","vX0"),location_est=c("FinalFound","ParamsInModel","vX0"),nameforplot="X0_BML",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,3),xlab="",ylab="",main="X(0)")
,list(param="B",evolmodel="mvslouch",location_true=c("modelparams","B"),location_est=c("FinalFound","ParamsInModel","B"),nameforplot="B_BML",f_cf=fL2_used,paramtype="",ylim=c(0,1),xlim=c(0,10),xlab="",ylab="",main="B")
,list(param="StS",evolmodel="bm",location_true=c("mvSLOUCHsummary","PointSummary","StS"),location_est=c("ParamSummary","StS"),nameforplot="StS",f_cf=fL2_used,paramtype="",ylim=c(0,1),xlim=c(0,5),xlab="",ylab="",main=expression(Sigma))
,list(param="StS",evolmodel="bm_indep",location_true=c("mvSLOUCHsummary","PointSummary","StS"),location_est=c("ParamSummary","StS"),nameforplot="StS",f_cf=fL2_used,paramtype="",ylim=c(0,1),xlim=c(0,5),xlab="",ylab="",main=expression(Sigma))
,list(param="StS",evolmodel="ouch",location_true=c("mvSLOUCHsummary","PointSummary","StS"),location_est=c("FinalFound","ParamSummary","StS"),nameforplot="StS",f_cf=fL2_used,paramtype="",ylim=c(0,1),xlim=c(0,5),xlab="",ylab="",main=expression(Sigma))
,list(param="StS",evolmodel="mvslouch",location_true=c("mvSLOUCHsummary","PointSummary","StS"),location_est=c("FinalFound","ParamSummary","StS"),nameforplot="StS_BML",f_cf=fL2_used,paramtype="",ylim=c(0,1),xlim=c(0,5),xlab="",ylab="",main=expression(Sigma))
,list(param="halflives",evolmodel="ouch",location_true=c("mvSLOUCHsummary","PointSummary","phyl.halflife","halflives"),location_est=c("FinalFound","ParamSummary","phyl.halflife","halflives"),nameforplot="halflives",f_cf=function(x,y){x<-x[2,];y<-y[2,];fL2(Re(x),Re(y))},paramtype="",ylim=c(0,0.5),xlim=c(0,1),xlab="",ylab="",main="half-life")
,list(param="halflives",evolmodel="mvslouch",location_true=c("mvSLOUCHsummary","PointSummary","phyl.halflife","halflives"),location_est=c("FinalFound","ParamSummary","phyl.halflife","halflives"),nameforplot="halflives_BML",f_cf=function(x,y){x<-x[2,];y<-y[2,];fL2(Re(x),Re(y))},paramtype="",ylim=c(0,0.5),xlim=c(0,50),xlab="",ylab="",main="half-life")
,list(param="eigenvalues",evolmodel="ouch",location_true=c("mvSLOUCHsummary","PointSummary","phyl.halflife","halflives"),location_est=c("FinalFound","ParamSummary","phyl.halflife","halflives"),nameforplot="eigenvalues",f_cf=function(x,y){x<-x[1,];y<-y[1,];fL2(Re(x),Re(y))},paramtype="",ylim=c(0,1),xlim=c(0,5),xlab="",ylab="",main="eigenvalue")
,list(param="eigenvalues",evolmodel="mvslouch",location_true=c("mvSLOUCHsummary","PointSummary","phyl.halflife","halflives"),location_est=c("FinalFound","ParamSummary","phyl.halflife","halflives"),nameforplot="eigenvalues_BML",f_cf=function(x,y){x<-x[1,];y<-y[1,];fL2(Re(x),Re(y))},paramtype="",ylim=c(0,1),xlim=c(0,5),xlab="",ylab="",main="eigenvalue")
,list(param="corr.matrix",evolmodel="ouch",location_true=c("mvSLOUCHsummary","PointSummary","corr.matrix"),location_est=c("FinalFound","ParamSummary","corr.matrix"),nameforplot="corrmatrix",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,1),xlab="",ylab="",main="correlation matrix")
,list(param="corr.matrix",evolmodel="mvslouch",location_true=c("mvSLOUCHsummary","PointSummary","corr.matrix"),location_est=c("FinalFound","ParamSummary","corr.matrix"),nameforplot="corrmatrix_BML",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,1),xlab="",ylab="",main="correlation matrix")
,list(param="stationary.corr.matrix",evolmodel="ouch",location_true=c("mvSLOUCHsummary","PointSummary","stationary.corr.matrix"),location_est=c("FinalFound","ParamSummary","stationary.corr.matrix"),nameforplot="statcorrmatrix",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,1),xlab="",ylab="",main="stationary correlation matrix")
,list(param="stationary.corr.matrix",evolmodel="mvslouch",location_true=c("mvSLOUCHsummary","PointSummary","stationary.corr.matrix"),location_est=c("FinalFound","ParamSummary","stationary.corr.matrix"),nameforplot="statcorrmatrix_BML",f_cf=fL2_used,paramtype="",ylim=c(0,2),xlim=c(0,1),xlab="",ylab="",main="stationary correlation matrix")
,list(param="trait.regression",evolmodel="ouch",location_true=c("mvSLOUCHsummary","PointSummary","trait.regression"),location_est=c("FinalFound","ParamSummary","trait.regression"),nameforplot="traitreg",f_cf=fL2list_used,paramtype="",ylim=c(0,2),xlim=c(0,0.1),xlab="",ylab="",main="trait regression")
,list(param="limiting.trait.regression",evolmodel="ouch",location_true=c("mvSLOUCHsummary","PointSummary","limiting.trait.regression"),location_est=c("FinalFound","ParamSummary","limiting.trait.regression"),nameforplot="limtraitreg",f_cf=fL2list_used,paramtype="",ylim=c(0,0.5),xlim=c(0,1),xlab="",ylab="",main="limiting trait regression")
,list(param="evolutionary.regression",evolmodel="mvslouch",location_true=c("mvSLOUCHsummary","PointSummary","evolutionary.regression"),location_est=c("FinalFound","ParamSummary","evolutionary.regression"),nameforplot="evolreg_BML",f_cf=fL2_used,paramtype="",ylim=c(0,0.5),xlim=c(0,1),xlab="",ylab="",main="evolutionary regression")
,list(param="optimal.regression",evolmodel="mvslouch",location_true=c("mvSLOUCHsummary","PointSummary","optimal.regression"),location_est=c("FinalFound","ParamSummary","optimal.regression"),nameforplot="optreg_BML",f_cf=fL2_used,paramtype="",ylim=c(0,0.5),xlim=c(0,1),xlab="",ylab="",main="optimal regression")
)


## create the boxplots 
lparamestsdists<-vector("list",length(lsimmodellist))

## create some dummy date in order to obtain parameter summaries
lsimmodellist<-sapply(1:length(lsimmodellist),function(i,lsimmodellist,dummy_data,dummy_phyltree,num_traits){
    n_dummy<-20
    k_dummy<-num_traits[i]
    tree_height<-1
    dummy_data<-matrix(rnorm(k_dummy*n_dummy),ncol=k_dummy,nrow=n_dummy)
    colnames(dummy_data)<-paste0("tr",1:k_dummy)
    dummy_phyltree<-TreeSim::sim.bd.taxa(n_dummy, 1, lambda=1, mu=0)[[1]]
    rownames(dummy_data)<-dummy_phyltree$tip.label ## so that mvSLOUCH summary functions do not raise warnings
    ## rescale to height 1 as in the simulations
    dummy_phyltree$edge.length<-tree_height*dummy_phyltree$edge.length/max(ape::node.depth.edgelength(dummy_phyltree))
    
    x<-lsimmodellist[[i]]
    ## phylogenies are all of height 1
    x$mvSLOUCHsummary<-switch(x$lmodelcf$evolmodel,
	bm=mvSLOUCH::SummarizeBM(phyltree=dummy_phyltree, mData=dummy_data, modelParams=x$modelparams, t = c(1)),
	bm_indep=mvSLOUCH::SummarizeBM(phyltree=dummy_phyltree, mData=dummy_data, modelParams=x$modelparams, t = c(1),dof=2*k_dummy),
	ouch=mvSLOUCH::SummarizeOUCH(phyltree=dummy_phyltree, mData=dummy_data, modelParams=x$modelparams, t = c(1),Atype=x$lmodelcf$Atype,Syytype=x$lmodelcf$Syytype),
	mvslouch=mvSLOUCH::SummarizeMVSLOUCH(phyltree=dummy_phyltree, mData=dummy_data, modelParams=x$modelparams, t = c(1),Atype=x$lmodelcf$Atype,Syytype=x$lmodelcf$Syytype)
    )[[1]] ## there are multiple time points
    x
},lsimmodellist=lsimmodellist,dummy_data=dummy_data,dummy_phyltree=dummy_phyltree,num_traits,simplify=FALSE)
## =================================================================================================================================



for (k in 1:length(lsimmodellist)){
    lparamestsdists[[k]]<-vector("list",length(lparamstocheck))
    names(lparamestsdists)[k]<-lsimmodellist[[k]]$modelnameforplot
    for (i in 1:length(lparamstocheck)){
	lparamestsdists[[k]][[i]]<-vector("list",length(vN))
	names(lparamestsdists[[k]][[i]])<-paste0("N_",vN)
	names(lparamestsdists[[k]])[i]<-paste0(lparamstocheck[[i]]$evolmodel,"_",lparamstocheck[[i]]$param)
	for (j in 1:length(vN)){
	    lparamestsdists[[k]][[i]][[j]]<-vector("list",1)
	    names(lparamestsdists[[k]][[i]][[j]])<-c("vdists")
	}
    }
}

dir.create(c_boxplotdir,"/", showWarnings = FALSE)

for(lmodel in lsimmodellist){
    dir.create(paste0(c_boxplotdir,"/",lmodel$modelnameforplot,"/"), showWarnings = FALSE)
    vNtmp<-vN
    if (lmodel$id=="05"){vNtmp<-c(32,64,128)}
    for (n in vNtmp){
	cfilename_lres<- paste0(c_dirpreffix,c_simuldir,"/",lmodel$id,"_ModelSelection/",c_fileprefix,"_",lmodel$id,"_N_",n,c_file_suffix,".RData")
	lparests<-f_getres(cfilename_lres,list(lmodel$lmodelcf),TRUE)
	for (paramtocheck in lparamstocheck){

	    b_paraminmodel<-paramtocheck$evolmodel==lmodel$lmodelcf$evolmodel	    
	    if (b_paraminmodel){
	    
		trueval<-f_findvalueinlist(lmodel,paramtocheck$location_true)
		vdists<-sapply(lparests,function(lpar,trueval,vlocation,f_cf){
		    estparamval<-f_findvalueinlist(lpar[[1]]$result,vlocation)
		    res<-NA
		    if (!is.null(estparamval)){
			res<-f_cf(trueval,estparamval)
		    }
		    res
		},vlocation=paramtocheck$location_est,trueval=trueval,f_cf=paramtocheck$f_cf,simplify=TRUE)
		k_modelindex<-which(names(lparamestsdists)==lmodel$modelnameforplot)
		i_paramindex<-which(names(lparamestsdists[[k_modelindex]])==paste0(paramtocheck$evolmodel,"_",paramtocheck$param))
		j_nindex<-which(names(lparamestsdists[[k_modelindex]][[i_paramindex]])==paste0("N_",n))
		lparamestsdists[[k_modelindex]][[i_paramindex]][[j_nindex]]$vdists<-vdists##vdists_orig

	    	
	    }	    	    
	}
    }

    
    k_modelindex<-which(names(lparamestsdists)==lmodel$modelnameforplot)
    for (paramtocheck in lparamstocheck){
        b_paraminmodel<-paramtocheck$evolmodel==lmodel$lmodelcf$evolmodel	    
        if (b_paraminmodel){
    	    i_paramindex<-which(names(lparamestsdists[[k_modelindex]])==paste0(paramtocheck$evolmodel,"_",paramtocheck$param))
    	    mdists<-c(NA,NA)
    	    for (n in vN){
    		j_nindex<-which(names(lparamestsdists[[k_modelindex]][[i_paramindex]])==paste0("N_",n))
    		mdists<-rbind(mdists,cbind(lparamestsdists[[k_modelindex]][[i_paramindex]][[j_nindex]]$vdists,rep(n,length(lparamestsdists[[k_modelindex]][[i_paramindex]][[j_nindex]]$vdists))))
    	    }
    	    mdists<-mdists[-1,,drop=FALSE]
    	    colnames(mdists)<-c("dist","n")
	    mdists<-as.data.frame(mdists)
	    sum_sq_dists<-sum(mdists$dist^2,na.rm=TRUE)
    	    if ((!is.na(sum_sq_dists))&&((sum_sq_dists)>0)){
		mdists$n<-as.factor(mdists$n)
		mdists$logdist<-mdists$dist 
		mdists$logdist[which(sapply(mdists$logdist,function(x){isTRUE(all.equal(x,0))},simplify=TRUE))]<-NA
		mdists$logdist<-log(mdists$logdist)	    
		vboxplotpdffilename<-paste0(c_boxplotdir,"/",lmodel$modelnameforplot,"/boxplot_",paramtocheck$nameforplot,"_",lmodel$modelnameforplot,".pdf")
		pdf(vboxplotpdffilename)
		plot(logdist~n,data=mdists,xlab="n",ylab=paramtocheck$ylab,main=paramtocheck$main,pch=19)
		dev.off()
	    }
	}
    }
}

