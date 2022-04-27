## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Rotation Invariance in non–Brownian motion Phylogenetic Comparative Methods: A comment on Adams and Collyer (2018)"
## It generates the values of Tab. 3 and should be run in the directory
## KJVJMRKK_mvSLOUCH/SimulationStudy/SimulationRuns/01_ModelSelection/
## containg the SimulationReestimation_SetupID_01_N_*.RData files


## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


vn<-c(32,64,128,256,512,1024,2048)

lalpha<-vector("list",length(vn))
for (r in 1:length(vn)){
    n<-vn[r]
    load(paste0("SimulationReestimation_SetupID_01_N_",n,".RData"))
    valpha<-rep(NA,length(lres))
    for (i in 1:length(lres)){
	BMAICc<-NA
	bestOUBM1AICc<-Inf
	jvalmvslouch<-NA
	for (j in 1:length(lres[[i]]$mvsl_estim_res$model.setups)){
	    if (lres[[i]]$mvsl_estim_res$model.setups[[j]]$evolmodel=="bm"){
		BMAICc<-lres[[i]]$mvsl_estim_res$testedModels[[j]]$aic.c
	    }else{
		if ((lres[[i]]$mvsl_estim_res$model.setups[[j]]$evolmodel=="mvslouch")&&(lres[[i]]$mvsl_estim_res$model.setups[[j]]$Atype=="SingleValueDiagonal")){		
		    tmpAICc<-lres[[i]]$mvsl_estim_res$testedModels[[j]]$aic.c
		    if ((!is.null(tmpAICc)) &&(!is.na(tmpAICc)) && (tmpAICc<bestOUBM1AICc)){
			bestOUBM1AICc<-tmpAICc;jvalmvslouch<-j
		    }
		}
	    }	    
	}
	if ((!is.null(BMAICc)) &&(!is.na(BMAICc)) && (!is.null(bestOUBM1AICc)) &&(!is.na(bestOUBM1AICc)) && (!is.na(bestOUBM1AICc)) && (bestOUBM1AICc<BMAICc)){
	    if (lres[[i]]$mvsl_estim_res$testedModels[[j]]$result$MaxLikFound[1]!="Same as final found"){
		valpha[i]<-lres[[i]]$mvsl_estim_res$testedModels[[j]]$result$MaxLikFound$ParamsInModel$A[1,1]
	    }else{
		valpha[i]<-lres[[i]]$mvsl_estim_res$testedModels[[j]]$result$FinalFound$ParamsInModel$A[1,1]
	    }	    
	}
    }
    lalpha[[r]]<-valpha
}

lvn<-length(vn)
print("average t_0.5")
to_print<-""
for(i in 1:lvn){to_print<-paste0(to_print,mean(log(2)/(2*lalpha[[i]]),na.rm=TRUE),", ")}
cat(paste0(to_print,"\n"))
print("==================================================")
print("variance t_0.5")
to_print<-""
for(i in 1:lvn){to_print<-paste0(to_print,var(log(2)/(2*lalpha[[i]]),na.rm=TRUE),", ")}
cat(paste0(to_print,"\n"))
print("==================================================")
print("% t_0.5 < 0")
to_print<-""
cutoffval<-0;for(i in 1:lvn){to_print<-paste0(to_print,100*sum(log(2)/(2*lalpha[[i]])<cutoffval,na.rm=TRUE)/sum(!is.na(lalpha[[i]]),na.rm=TRUE),", ")}
cat(paste0(to_print,"\n"))
print("==================================================")
print("% t_0.5 > 0.5")
to_print<-""
cutoffval<-0.5;for(i in 1:lvn){to_print<-paste0(to_print,100*sum(log(2)/(2*lalpha[[i]])>cutoffval,na.rm=TRUE)/sum(!is.na(lalpha[[i]]),na.rm=TRUE),", ")}
cat(paste0(to_print,"\n"))
print("==================================================")
print("% t_0.5 > 1")
to_print<-""
cutoffval<-1;for(i in 1:lvn){to_print<-paste0(to_print,100*sum(log(2)/(2*lalpha[[i]])>cutoffval,na.rm=TRUE)/sum(!is.na(lalpha[[i]]),na.rm=TRUE),", ")}
cat(paste0(to_print,"\n"))
print("==================================================")
