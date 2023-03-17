## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Rotation Invariance in non–Brownian motion Phylogenetic Comparative Methods: A comment on Adams and Collyer (2018)"
## It generates the values of Tab. 4,5 and should be run in the main directory KJVJMRKK_mvSLOUCH/SimulationStudy/
## which contains the SimulationRuns/ directory


## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .






f_getalphavalues<-function(vn,sim_model,cf_model,setupid,dir_prefix=""){
    lalpha<-vector("list",length(vn))
    for (r in 1:length(vn)){
	n<-vn[r]
        load(paste0(dir_prefix,"SimulationReestimation_SetupID_",setupid,"_N_",n,".RData"))
	valpha<-rep(NA,length(lres))
        for (i in 1:length(lres)){
	    BMAICc<-NA
	    bestOUBM1AICc<-Inf
	    jvalmvslouch<-NA
	    for (j in 1:length(lres[[i]]$mvsl_estim_res$model.setups)){
		if (lres[[i]]$mvsl_estim_res$model.setups[[j]]$evolmodel==sim_model$evolmodel){ 
		    BMAICc<-lres[[i]]$mvsl_estim_res$testedModels[[j]]$aic.c
		}else{
		    if ((!is.null(lres[[i]]$mvsl_estim_res$model.setups[[j]]$evolmodel))&&(lres[[i]]$mvsl_estim_res$model.setups[[j]]$evolmodel==cf_model$evolmodel)&&(!is.null(lres[[i]]$mvsl_estim_res$model.setups[[j]]$Atype))&&(lres[[i]]$mvsl_estim_res$model.setups[[j]]$Atype==cf_model$Atype)&&
		    ( (is.null(cf_model$diagA)&&is.null(lres[[i]]$mvsl_estim_res$model.setups[[j]]$diagA))
		    || ((!is.null(cf_model$diagA))&&(!is.null(lres[[i]]$mvsl_estim_res$model.setups[[j]]$diagA))&&(lres[[i]]$mvsl_estim_res$model.setups[[j]]$diagA==cf_model$diagA)
		    )
		    )){ 
			tmpAICc<-lres[[i]]$mvsl_estim_res$testedModels[[j]]$aic.c
			if ((!is.null(tmpAICc)) &&(!is.na(tmpAICc)) && (tmpAICc<bestOUBM1AICc)){
			    bestOUBM1AICc<-tmpAICc;jvalmvslouch<-j
			}
		    }		    
		}	    
	    }
	    if ((!is.null(BMAICc)) &&(!is.na(BMAICc)) && (!is.null(bestOUBM1AICc)) &&(!is.na(bestOUBM1AICc)) && (!is.na(bestOUBM1AICc)) && (bestOUBM1AICc<BMAICc)){
		if (lres[[i]]$mvsl_estim_res$testedModels[[jvalmvslouch]]$result$MaxLikFound[1]!="Same as final found"){
		    valpha[i]<-lres[[i]]$mvsl_estim_res$testedModels[[jvalmvslouch]]$result$MaxLikFound$ParamsInModel$A[1,1]
		}else{
		    valpha[i]<-lres[[i]]$mvsl_estim_res$testedModels[[jvalmvslouch]]$result$FinalFound$ParamsInModel$A[1,1]
		}	    
	    }
	}
	lalpha[[r]]<-valpha
    }
    lalpha
}


f_cfhalflives<-function(vn,sim_model,cf_model,setupid,dir_prefix=""){
    lalpha<-f_getalphavalues(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)
    lvn<-length(vn)
    print("Model simulated under:")
    print(sim_model)
    print("Model being compared with:")
    print(cf_model)
    print("average t_0.5")
    to_print<-""
    for(i in 1:lvn){to_print<-paste0(to_print,round(mean(log(2)/(lalpha[[i]]),na.rm=TRUE),3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("variance t_0.5")
    to_print<-""
    for(i in 1:lvn){to_print<-paste0(to_print,round(var(log(2)/(lalpha[[i]]),na.rm=TRUE),3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("median t_0.5")
    to_print<-""
    for(i in 1:lvn){to_print<-paste0(to_print,round(median(log(2)/(lalpha[[i]]),na.rm=TRUE),3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("lower quartile t_0.5")
    to_print<-""
    for(i in 1:lvn){to_print<-paste0(to_print,round(quantile(log(2)/(lalpha[[i]]),na.rm=TRUE)[2],3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("upper quartile t_0.5")
    to_print<-""
    for(i in 1:lvn){to_print<-paste0(to_print,round(quantile(log(2)/(lalpha[[i]]),na.rm=TRUE)[4],3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("% t_0.5 < 0")
    to_print<-""
    cutoffval<-0;for(i in 1:lvn){to_print<-paste0(to_print,round(100*sum(log(2)/(lalpha[[i]])<cutoffval,na.rm=TRUE)/sum(!is.na(lalpha[[i]]),na.rm=TRUE),3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("% t_0.5 > 0.5")
    to_print<-""
    cutoffval<-0.5;for(i in 1:lvn){to_print<-paste0(to_print,round(100*sum(log(2)/(lalpha[[i]])>cutoffval,na.rm=TRUE)/sum(!is.na(lalpha[[i]]),na.rm=TRUE),3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("% t_0.5 > 1")
    to_print<-""
    cutoffval<-1;for(i in 1:lvn){to_print<-paste0(to_print,round(100*sum(log(2)/(lalpha[[i]])>cutoffval,na.rm=TRUE)/sum(!is.na(lalpha[[i]]),na.rm=TRUE),3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("% t_0.5 > 2")
    to_print<-""
    cutoffval<-2;for(i in 1:lvn){to_print<-paste0(to_print,round(100*sum(log(2)/(lalpha[[i]])>cutoffval,na.rm=TRUE)/sum(!is.na(lalpha[[i]]),na.rm=TRUE),3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
    print("% t_0.5 > 3")
    to_print<-""
    cutoffval<-3;for(i in 1:lvn){to_print<-paste0(to_print,round(100*sum(log(2)/(lalpha[[i]])>cutoffval,na.rm=TRUE)/sum(!is.na(lalpha[[i]]),na.rm=TRUE),3),", ")}
    cat(paste0(to_print,"\n"))
    print("==================================================")
}

sink("cf_halflives_bmindep.txt")

sim_model<-list(evolmodel="bm_indep")

print("Setup 01 p=4")
dir_prefix<-"SimulationRuns/01_ModelSelection/"
setupid="01"
vn<-c(32,64,128,256,512,1024,2048)

cf_model<-list(evolmodel="ouch","Atype"="SingleValueDiagonal","diagA"="Positive")
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)

cf_model<-list(evolmodel="ouch","Atype"="SingleValueDiagonal","diagA"=NULL)
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

print("Setup 02 p=8")
dir_prefix<-"SimulationRuns/02_ModelSelection/"
setupid="02"
vn<-c(32,64,128,256,512,1024,2048)

cf_model<-list(evolmodel="ouch","Atype"="SingleValueDiagonal","diagA"="Positive")
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)

cf_model<-list(evolmodel="ouch","Atype"="SingleValueDiagonal","diagA"=NULL)
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

print("Setup 03 p=12")
dir_prefix<-"SimulationRuns/03_ModelSelection/"
setupid="03"
vn<-c(32,64,128,256,512,1024,2048)

cf_model<-list(evolmodel="ouch","Atype"="SingleValueDiagonal","diagA"="Positive")
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)

cf_model<-list(evolmodel="ouch","Atype"="SingleValueDiagonal","diagA"=NULL)
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

print("Setup 04 p=16")
dir_prefix<-"SimulationRuns/04_ModelSelection/"
setupid="04"
vn<-c(32,64,128,256,512,1024,2048)

cf_model<-list(evolmodel="ouch","Atype"="SingleValueDiagonal","diagA"="Positive")
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)

cf_model<-list(evolmodel="ouch","Atype"="SingleValueDiagonal","diagA"=NULL)
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

print("Setup 05 p=32")
dir_prefix<-"SimulationRuns/05_ModelSelection/"
#vn<-c(32,64,128,256,512,1024,2048)
setupid="05"
vn<-c(32)#,64,128)

cf_model<-list(evolmodel="ouch","Atype"="DecomposablePositive","diagA"="Positive")
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)

cf_model<-list(evolmodel="ouch","Atype"="DecomposablePositive","diagA"=NULL)
f_cfhalflives(vn=vn,sim_model=sim_model,cf_model=cf_model,setupid=setupid,dir_prefix=dir_prefix)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

sink()
