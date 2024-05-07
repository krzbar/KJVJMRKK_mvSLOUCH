## This file accompanies the manuscript: 
## Bartoszek, Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## The code here generates the output for 
## Appendix SG: mvSLOUCH timings

library(parallel)
source("tosource_for_mvSLOUCH.R")


numcores<-8
b_use_random_seed_from_manuscript<-FALSE
b_should_random_seed_be_saved<-TRUE
filename_suffix<-"" ## initial _ has to be here
filename_prefix<-"mvSLOUCH_newold_timings"
results_directory<-"./mvSLOUCH_timings/"
tmp_directory<-paste0(results_directory,"tmp_res/")
BMtextresfile<-paste0(results_directory,filename_prefix,"_BM",filename_suffix,".txt")
OUOUtextresfile<-paste0(results_directory,filename_prefix,"_OUOU",filename_suffix,".txt")
OUBMtextresfile<-paste0(results_directory,filename_prefix,"_OUBM",filename_suffix,".txt")
outfiletimingresults<-paste0(results_directory,filename_prefix,"_summary",filename_suffix,".txt")

num_repeats<-32
vN<-c(5,10,25,50,100,200,400,500)
BMparams<-list(vX0=matrix(0,nrow=4,ncol=1),Sxx=rbind(c(1,0.25,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)))
OUOUparams<-list(vY0=matrix(c(1,-1,0,0),ncol=1,nrow=4),A=rbind(c(9,3,2,-2),c(2,5,-2,2),c(0,0,5,1),c(0,0,3,2)),mPsi=matrix(c(0,0,0,0),ncol=1,nrow=4),Syy=rbind(c(1,0.25,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)))
OUBMparams<-list(vY0=matrix(c(1,-1),ncol=1,nrow=2),A=rbind(c(9,3),c(2,5)),B=cbind(c(2,-2),c(-2,2)),mPsi=matrix(c(0,0),ncol=1,nrow=2),Syy=rbind(c(1,0.25),c(0,1)),vX0=matrix(c(0,0),ncol=1,nrow=2),Sxx=diag(c(1,1)),Syx=matrix(0,ncol=2,nrow=2),Sxy=matrix(0,ncol=2,nrow=2))



dir.create(results_directory,showWarnings=FALSE)
dir.create(tmp_directory,showWarnings=FALSE)
if (!exists("cl")){cl <- makeCluster(getOption("cl.cores", numcores),outfile="")}
f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,results_directory,filename_prefix,"RandomSeed_Start",filename_suffix)
lresults_mvSLOUCH_timings<-sapply(vN,function(n,results_directory,filename_prefix,filename_suffix,b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,num_repeats,BMparams,OUOUparams,OUBMparams,BMtextresfile,OUOUtextresfile,OUBMtextresfile,cl,outfiletimingresults){
    f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,results_directory,filename_prefix,paste0("RandomSeed_Start_n_",n),filename_suffix)
    lresults_mvSLOUCH_timings_n<-parSapply(cl,1:num_repeats,function(i,n,results_directory,filename_prefix,filename_suffix,b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,num_repeats,BMparams,OUOUparams,OUBMparams,BMtextresfile,OUOUtextresfile,OUBMtextresfile){
	source("tosource_for_mvSLOUCH.R")
	
        lres_n_i<-vector("list",6)
        names(lres_n_i)<-c("timings","ape_phyltree","ouch_phyltree","BM","OUOU","OUBM")
        v_timings<-rep(NA,6)
        names(v_timings)<-c("BM_old","BM_new","OUOU_old","OUOU_new","OUBM_old","OUBM_new")
        lres_n_i$BM<-vector("list",7)
        names(lres_n_i$BM)<-c("timings_se","ape_phyltree","ouch_phyltree","data","BMold","BMnew","BMparams")
	lres_n_i$OUOU<-vector("list",7)
        names(lres_n_i$OUOU)<-c("timings_se","ape_phyltree","ouch_phyltree","data","OUOUold","OUOUnew","OUOUparams")
	lres_n_i$OUBM<-vector("list",7)
        names(lres_n_i$OUBM)<-c("timings_se","ape_phyltree","ouch_phyltree","data","OUBMold","OUBMnew","OUBMparams")

    
## Simulate all of the data
        f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,filename_prefix,paste0("RandomSeed_Start_n_",n,"_iter_",i,"_start"),filename_suffix)
	ape_phyltree<-TreeSim::sim.bd.taxa(n, 1, lambda=1, mu=0)[[1]]
        ## rescale tree to height 1 so all simulations and estimations are comparable
    	ape_phyltree$edge.length<-1*ape_phyltree$edge.length/max(ape::node.depth.edgelength(ape_phyltree))
	
	lres_n_i$ape_phyltree<-ape_phyltree
        
	BMdata<-mvSLOUCH::simulBMProcPhylTree(ape_phyltree,X0=BMparams$vX0,Sigma=BMparams$Sxx)
        BMdata<-BMdata[ape_phyltree$tip.label,,drop=FALSE]
        lres_n_i$BM$data<-BMdata

	OUOUdata<-mvSLOUCH::simulOUCHProcPhylTree(ape_phyltree,OUOUparams)
        OUOUdata<-OUOUdata[ape_phyltree$tip.label,,drop=FALSE]
        lres_n_i$OUOU$data<-OUOUdata

	OUBMdata<-mvSLOUCH::simulMVSLOUCHProcPhylTree(ape_phyltree,OUBMparams)
        OUBMdata<-OUBMdata[ape_phyltree$tip.label,,drop=FALSE]
        lres_n_i$OUBM$data<-OUBMdata
## =================================================================================================================

## BM
        lres_n_i$BM$timings_se<-matrix(NA,2,2)
        colnames(lres_n_i$BM$timings_se)<-c("start","end")
        rownames(lres_n_i$BM$timings_se)<-c("old","new")
        lres_n_i$BM$BMparams<-BMparams
        lres_n_i$BM$ape_phyltree<-ape_phyltree

        f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,filename_prefix,paste0("RandomSeed_Start_n_",n,"_iter_",i,"_oldBM"),filename_suffix)
	stime<-as.POSIXct(NA);etime<-as.POSIXct(NA);BMoldres<-NA;ouch_phyltree<-NA
	tryCatch({
	    stime<-Sys.time()
	    ouch_phyltree<-ouch::ape2ouch(ape_phyltree,scale=1)     
    	    ### Correct the names of the internal node labels
    	    ouch_phyltree@nodelabels[1:(ouch_phyltree@nnodes-ouch_phyltree@nterm)]<-as.character(1:(ouch_phyltree@nnodes-ouch_phyltree@nterm))
	    BMoldres<-mvSLOUCHold::BrownianMotionModel(ouch_phyltree,BMdata)
	    etime<-Sys.time()
	},error=function(e){})
	lres_n_i$BM$timings_se[1,]<-c(stime,etime)
        lres_n_i$BM$ouch_phyltree<-ouch_phyltree
        lres_n_i$BM$BMold<-BMoldres
	v_timings[1]<-as.numeric(etime-stime,units="secs")
	sink(BMtextresfile,append=TRUE)
	cat(paste0("BM n ",n," i ",i," old: ",v_timings[1],"\n"))
	sink()

	f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,filename_prefix,paste0("RandomSeed_Start_n_",n,"_iter_",i,"_newBM"),filename_suffix)
	stime<-as.POSIXct(NA);etime<-as.POSIXct(NA);BMnewres<-NA;ouch_phyltree<-NA
	tryCatch({
	    stime<-Sys.time()
    	    BMnewres<-mvSLOUCH::BrownianMotionModel(ape_phyltree,BMdata)
	    etime<-Sys.time()
	},error=function(e){})
	lres_n_i$BM$timings_se[1,]<-c(stime,etime)
        lres_n_i$BM$BMnew<-BMnewres
	v_timings[2]<-as.numeric(etime-stime,units="secs")
	sink(BMtextresfile,append=TRUE)
	cat(paste0("BM n ",n," i ",i," new: ",v_timings[2],"\n"))
	sink()
        BMoldnew<-lres_n_i$BM
	save(BMoldnew,file=paste0(tmp_directory,filename_prefix,"BMres_n_",n,"_iter_",i,filename_suffix,".RData"))
## =================================================================================================================


## OUOU
        lres_n_i$OUOU$timings_se<-matrix(NA,2,2)
        colnames(lres_n_i$OUOU$timings_se)<-c("start","end")
        rownames(lres_n_i$OUOU$timings_se)<-c("old","new")
        lres_n_i$OUOU$OUOUparams<-OUOUparams
        lres_n_i$OUOU$ape_phyltree<-ape_phyltree

        f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,filename_prefix,paste0("RandomSeed_Start_n_",n,"_iter_",i,"_oldOUOU"),filename_suffix)
	stime<-as.POSIXct(NA);etime<-as.POSIXct(NA);OUOUoldres<-NA;ouch_phyltree<-NA
	tryCatch({
	    stime<-Sys.time()
	    ouch_phyltree<-ouch::ape2ouch(ape_phyltree,scale=1)     
    	    ### Correct the names of the internal node labels
    	    ouch_phyltree@nodelabels[1:(ouch_phyltree@nnodes-ouch_phyltree@nterm)]<-as.character(1:(ouch_phyltree@nnodes-ouch_phyltree@nterm))
	    OUOUoldres<-mvSLOUCHold::ouchModel(ouch_phyltree,OUOUdata,NULL,Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive")
	    etime<-Sys.time()
	},error=function(e){})
	lres_n_i$OUOU$timings_se[1,]<-c(stime,etime)
        lres_n_i$OUOU$ouch_phyltree<-ouch_phyltree
        lres_n_i$OUOU$OUOUold<-OUOUoldres
	v_timings[3]<-as.numeric(etime-stime,units="secs")
	sink(OUOUtextresfile,append=TRUE)
	cat(paste0("OUOU n ",n," i ",i," old: ",v_timings[3],"\n"))
	sink()

	f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,filename_prefix,paste0("RandomSeed_Start_n_",n,"_iter_",i,"_newOUOU"),filename_suffix)
	stime<-as.POSIXct(NA);etime<-as.POSIXct(NA);OUOUnewres<-NA;ouch_phyltree<-NA
	tryCatch({
	    stime<-Sys.time()
    	    OUOUnewres<-mvSLOUCH::ouchModel(ape_phyltree,OUOUdata,NULL,Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA="Positive")
	    etime<-Sys.time()
	},error=function(e){})
	lres_n_i$OUOU$timings_se[1,]<-c(stime,etime)
        lres_n_i$OUOU$OUOUnew<-OUOUnewres
	v_timings[4]<-as.numeric(etime-stime,units="secs")
	sink(OUOUtextresfile,append=TRUE)
	cat(paste0("OUOU n ",n," i ",i," new: ",v_timings[4],"\n"))
	sink()
        OUOUoldnew<-lres_n_i$OUOU
	save(OUOUoldnew,file=paste0(tmp_directory,filename_prefix,"OUOUres_n_",n,"_iter_",i,filename_suffix,".RData"))

## =================================================================================================================

## OUBM
        lres_n_i$OUBM$timings_se<-matrix(NA,2,2)
        colnames(lres_n_i$OUBM$timings_se)<-c("start","end")
        rownames(lres_n_i$OUBM$timings_se)<-c("old","new")
        lres_n_i$OUBM$OUBMparams<-OUBMparams
        lres_n_i$OUBM$ape_phyltree<-ape_phyltree

        f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,filename_prefix,paste0("RandomSeed_Start_n_",n,"_iter_",i,"_oldOUBM"),filename_suffix)
	stime<-as.POSIXct(NA);etime<-as.POSIXct(NA);OUBMoldres<-NA;ouch_phyltree<-NA
	tryCatch({
	    stime<-Sys.time()
	    ouch_phyltree<-ouch::ape2ouch(ape_phyltree,scale=1)     
    	    ### Correct the names of the internal node labels
    	    ouch_phyltree@nodelabels[1:(ouch_phyltree@nnodes-ouch_phyltree@nterm)]<-as.character(1:(ouch_phyltree@nnodes-ouch_phyltree@nterm))
	    OUBMoldres<-mvSLOUCHold::mvslouchModel(ouch_phyltree,OUBMdata,2,NULL,Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive")
	    etime<-Sys.time()
	},error=function(e){})
	lres_n_i$OUBM$timings_se[1,]<-c(stime,etime)
        lres_n_i$OUBM$ouch_phyltree<-ouch_phyltree
        lres_n_i$OUBM$OUBMold<-OUBMoldres
	v_timings[5]<-as.numeric(etime-stime,units="secs")
	sink(OUBMtextresfile,append=TRUE)
	cat(paste0("OUBM n ",n," i ",i," old: ",v_timings[5],"\n"))
	sink()

	f_loadsave_randomseed(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,tmp_directory,filename_prefix,paste0("RandomSeed_Start_n_",n,"_iter_",i,"_newOUBM"),filename_suffix)
	stime<-as.POSIXct(NA);etime<-as.POSIXct(NA);OUBMnewres<-NA;ouch_phyltree<-NA
	tryCatch({
	    stime<-Sys.time()
    	    OUBMnewres<-mvSLOUCH::mvslouchModel(ape_phyltree,OUBMdata,2,NULL,Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive")
	    etime<-Sys.time()
	},error=function(e){})
	lres_n_i$OUBM$timings_se[1,]<-c(stime,etime)
        lres_n_i$OUBM$OUBMnew<-OUBMnewres
	v_timings[6]<-as.numeric(etime-stime,units="secs")
	sink(OUBMtextresfile,append=TRUE)
	cat(paste0("OUBM n ",n," i ",i," new: ",v_timings[6],"\n"))
	sink()
        OUBMoldnew<-lres_n_i$OUBM
        save(OUBMoldnew,file=paste0(tmp_directory,filename_prefix,"OUBMres_n_",n,"_iter_",i,filename_suffix,".RData"))
	
	
	lres_n_i$timings<-v_timings	
	lres_n_i
## =================================================================================================================
    },n=n,results_directory=results_directory,filename_prefix=filename_prefix,filename_suffix=filename_suffix,b_use_random_seed_from_manuscript=b_use_random_seed_from_manuscript,b_should_random_seed_be_saved=b_should_random_seed_be_saved,tmp_directory=tmp_directory,num_repeats=num_repeats,BMparams=BMparams,OUOUparams=OUOUparams,OUBMparams=OUBMparams,BMtextresfile=BMtextresfile,OUOUtextresfile=OUOUtextresfile,OUBMtextresfile=OUBMtextresfile,simplify=FALSE)

    mTimings<-sapply(lresults_mvSLOUCH_timings_n,function(x){x$timings},simplify=TRUE)
    vMeans<-NA;vMeans<-tryCatch(apply(mTimings,1,mean,na.rm=TRUE),error=function(e){})
    vVars<-NA;vVars<-tryCatch(apply(mTimings,1,var,na.rm=TRUE),error=function(e){})
    vNumNAs<-apply(mTimings,1,function(x){sum(is.na(x))})
    sink(outfiletimingresults,append=TRUE)
    cat(paste0("n = ",n,"\n"))
    cat("Mean times \n")
    print(vMeans);cat("\n")
    cat("Var times \n")
    print(vVars);cat("\n")
    cat("Number of failed estimations \n")
    print(vNumNAs) ;cat("\n")
    cat("=================================================== \n")
    sink()
    lresults_mvSLOUCH_timings_n$mTimings<-mTimings
    lresults_mvSLOUCH_timings_n$vMeans<-vMeans
    lresults_mvSLOUCH_timings_n$vVars<-vVars
    lresults_mvSLOUCH_timings_n$vNumNAs<-vNumNAs
    save(lresults_mvSLOUCH_timings_n,file=paste0(results_directory,filename_prefix,"mvSLOUCHnewoldres_n_",n,filename_suffix,".RData"))
},results_directory=results_directory,filename_prefix=filename_prefix,filename_suffix=filename_suffix,b_use_random_seed_from_manuscript=b_use_random_seed_from_manuscript,b_should_random_seed_be_saved=b_should_random_seed_be_saved,tmp_directory=tmp_directory,num_repeats=num_repeats,BMparams=BMparams,OUOUparams=OUOUparams,OUBMparams=OUBMparams,BMtextresfile=BMtextresfile,OUOUtextresfile=OUOUtextresfile,OUBMtextresfile=OUBMtextresfile,cl=cl,outfiletimingresults=outfiletimingresults,simplify=FALSE)
