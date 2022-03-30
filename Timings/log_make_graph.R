## This file accompanies the manuscript: 
## Bartoszek, Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## User controlled variables ===========================================================================
filename_suffix<-"" ## initial _ has to be here
filename_prefix<-"mvSLOUCH_newold_timings"
results_directory<-"./mvSLOUCH_timings/"
tmp_directory<-paste0(results_directory,"tmp_res/")
vN<-c(5,10,25,50,100,200,400,500) ## numbers of tips for plotting
vcolours<-c("black","red","blue","green","orange","brown") ## these are for BM, OUOU, OUBM with different mvSLOUCH versions, order as in legend
## end of user controlled variables
## =====================================================================================================

l_timings<-vector("list",length(vN))
names(l_timings)<-paste0("n_",vN)

for (n in vN){
    l_timings[[which(names(l_timings)==paste0("n_",n))]]<-vector("list",4)
    names(l_timings[[which(names(l_timings)==paste0("n_",n))]])<-c("mTimings","vMeans","vVars","vNumNAs")
    load(paste0(results_directory,filename_prefix,"mvSLOUCHnewoldres_n_",n,filename_suffix,".RData"))
    l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings<-(lresults_mvSLOUCH_timings_n$mTimings)
}

plot(NA,xlim=c(0,max(vN)+1),ylim=c(-5,(max(apply(log(l_timings[[which(names(l_timings)==paste0("n_",max(vN)))]]$mTimings[c(1,3,5),,drop=FALSE]),1,median,na.rm=TRUE),na.rm=TRUE))),main="",xlab="n",ylab="log(time[s])",cex.lab=1.25,cex.axis=1.25)

v_pred_n<-seq(5,max(vN)+1,by=0.1)
for(i in 1:6){
    vMedians<-rep(NA,length(vN))
    names(vMedians)<-paste0("n_",vN)
    for (n in vN){
    len_bar<-0.025
        median_time<-median(log(l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings[i,]),na.rm=TRUE)
        vMedians[which(names(vMedians)==paste0("n_",n))]<-median((l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings[i,]),na.rm=TRUE)
    interquant_time<-quantile(log(l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings[i,]),na.rm=TRUE)[c(2,4)]
    interquant_time[1]<-median_time-1.5*(median_time-interquant_time[1])
    interquant_time[2]<-median_time+1.5*(interquant_time[2]-median_time)
        points(n,median_time,pch=19,col=vcolours[i])
    }
    sink("LineFitCoeffs.txt",append=TRUE)
    print(paste0("Setup ",i))
    if (i%%2==1){## we have old
    res_quadlm<-lm(((vMedians))~vN+I(vN^2))
    coeff_a<-res_quadlm$coefficients[3] #3  #ax^2+b+c [1]->c, [2]->b, [3]->a
        coeff_b<-res_quadlm$coefficients[2]
    coeff_c<-res_quadlm$coefficients[1]
    print(c(coeff_a,coeff_b,coeff_c))
    v_predtime<-log((coeff_a*v_pred_n^2+coeff_b*v_pred_n+coeff_c))
    }else{## we have new
        res_quadlm<-lm(((vMedians))~vN)
        coeff_b<-res_quadlm$coefficients[2]
    coeff_c<-res_quadlm$coefficients[1]
    v_predtime<-log((coeff_b*v_pred_n+coeff_c))
    print(c(coeff_b,coeff_c))
    }
    print("=====================================================")
    sink()
    if (i==1){points(v_pred_n[which(v_pred_n>150)],v_predtime[which(v_pred_n>150)],type="l",cex=0.5,col=vcolours[i],lwd=1)}
    else{points(v_pred_n,v_predtime,type="l",cex=0.5,col=vcolours[i],lwd=1)}
}

legend("bottomright",legend=c("BM 1.3.4","BM 2.7.3","OUOU 1.3.4", "OUOU 2.7.3", "OUBM 1.3.4", "OUBM 2.7.3"),col=vcolours,pch=19,bty="n")
