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
vN<-c(5,10,25,50,100,200,400,500) ## number of tips for plotting 
vcolours<-c("gray","red","blue","green","orange","brown") ## these are for BM, OUOU, OUBM with different mvSLOUCH versions, order as in legend
v_inside_graph<-c(2,4,6) ## which curves are to be in the small inside graph, order as in legend
## end of user controlled variables
## =====================================================================================================
l_timings<-vector("list",length(vN))
names(l_timings)<-paste0("n_",vN)
v_outside_graph<-setdiff(1:6,v_inside_graph)

for (n in vN){
    l_timings[[which(names(l_timings)==paste0("n_",n))]]<-vector("list",4)
    names(l_timings[[which(names(l_timings)==paste0("n_",n))]])<-c("mTimings","vMeans","vVars","vNumNAs")
    load(paste0(results_directory,filename_prefix,"mvSLOUCHnewoldres_n_",n,filename_suffix,".RData"))
    l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings<-(lresults_mvSLOUCH_timings_n$mTimings)
}

plot(NA,xlim=c(0,max(vN)+1),ylim=c(0,(max(l_timings[[which(names(l_timings)==paste0("n_",max(vN)))]]$mTimings[v_outside_graph,],na.rm=TRUE))),main="",xlab="n",ylab="time[s]")

v_pred_n<-seq(5,max(vN)+1,by=0.1)
for(i in v_outside_graph){
    vMedians<-rep(NA,length(vN))
    names(vMedians)<-paste0("n_",vN)
    for (n in vN){
    len_bar<-0.025
    #https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
        median_time<-median((l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings[i,]),na.rm=TRUE)
        vMedians[which(names(vMedians)==paste0("n_",n))]<-median_time
    interquant_time<-quantile((l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings[i,]),na.rm=TRUE)[c(2,4)]
    interquant_time[1]<-median_time-1.5*(median_time-interquant_time[1])
    interquant_time[2]<-median_time+1.5*(interquant_time[2]-median_time)
    segments(n,interquant_time[1], n, interquant_time[2], col=vcolours[i],lwd=2)
    segments(n-len_bar,interquant_time[1], n+len_bar, interquant_time[1], col=vcolours[i],lwd=2)
    segments(n-len_bar,interquant_time[2], n+len_bar, interquant_time[2], col=vcolours[i],lwd=2)
    }
    sink("LineFitCoeffs.txt",append=TRUE)
    print(paste0("Setup ",i))
    if (i%%2==1){## we have old
        res_quadlm<-lm((vMedians)~vN+I(vN^2))
    coeff_a<-res_quadlm$coefficients[3]  #ax^2+b+c [1]->c, [2]->b, [3]->a
        coeff_b<-res_quadlm$coefficients[2]
    coeff_c<-res_quadlm$coefficients[1]
    print(c(coeff_a,coeff_b,coeff_c))
    v_predtime<-(coeff_a*v_pred_n^2+coeff_b*v_pred_n+coeff_c)
    }else{## we have new
        res_quadlm<-lm((vMedians)~vN)
        coeff_b<-res_quadlm$coefficients[2]
    coeff_c<-res_quadlm$coefficients[1]
    v_predtime<-(coeff_b*v_pred_n+coeff_c)	
    print(c(coeff_b,coeff_c))
    }
    print("=====================================================")
    sink()
    points(v_pred_n,v_predtime,type="l",cex=0.5,col=vcolours[i],lwd=1)
}

legend("topleft",legend=c("BM 1.3.4","BM 2.7.3","OUOU 1.3.4", "OUOU 2.7.3", "OUBM 1.3.4", "OUBM 2.7.3"),col=vcolours,pch=19,bty="n")
par(fig=c(0.325,0.9, 0.45,1), new=TRUE, las=1, ps=9,bty="n")
plot(NA,xlim=c(0,max(vN)+1),ylim=c(0,max(apply(l_timings[[which(names(l_timings)==paste0("n_",max(vN)))]]$mTimings[v_inside_graph,,drop=FALSE],1,median,na.rm=TRUE),na.rm=TRUE)+50),main="",xlab="n",ylab="time[s]")

v_pred_n<-seq(5,max(vN)+1,by=0.1)
for(i in v_inside_graph){
    vMedians<-rep(NA,length(vN))
    names(vMedians)<-paste0("n_",vN)
    for (n in vN){
    len_bar<-0.025
    #https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
        median_time<-median((l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings[i,]),na.rm=TRUE)
        vMedians[which(names(vMedians)==paste0("n_",n))]<-median_time
    interquant_time<-quantile((l_timings[[which(names(l_timings)==paste0("n_",n))]]$mTimings[i,]),na.rm=TRUE)[c(2,4)]
    interquant_time[1]<-median_time-1.5*(median_time-interquant_time[1])
    interquant_time[2]<-median_time+1.5*(interquant_time[2]-median_time)
        points(n,median_time,pch=19,col=vcolours[i])
    lwd_use<-2
    if (i==1){lwd_use<-3}
    segments(n,interquant_time[1], n, interquant_time[2], col=vcolours[i],lwd=lwd_use)
    segments(n-len_bar,interquant_time[1], n+len_bar, interquant_time[1], col=vcolours[i],lwd=lwd_use)
    segments(n-len_bar,interquant_time[2], n+len_bar, interquant_time[2], col=vcolours[i],lwd=lwd_use)
    }
    sink("LineFitCoeffs.txt",append=TRUE)
    print(paste0("Setup ",i))
    if (i%%2==1){## we have old
        res_quadlm<-lm((vMedians)~vN+I(vN^2))
    coeff_a<-res_quadlm$coefficients[3]  #ax^2+b+c [1]->c, [2]->b, [3]->a
        coeff_b<-res_quadlm$coefficients[2]
    coeff_c<-res_quadlm$coefficients[1]
    print(c(coeff_a,coeff_b,coeff_c))
    v_predtime<-(coeff_a*v_pred_n^2+coeff_b*v_pred_n+coeff_c)
    }else{## we have new
        res_quadlm<-lm((vMedians)~vN)
        coeff_b<-res_quadlm$coefficients[2]
    coeff_c<-res_quadlm$coefficients[1]
    v_predtime<-(coeff_b*v_pred_n+coeff_c)	
    print(c(coeff_b,coeff_c))
    }
    print("=====================================================")
    sink()
    points(v_pred_n,v_predtime,type="l",cex=0.5,col=vcolours[i],lwd=1)
}
