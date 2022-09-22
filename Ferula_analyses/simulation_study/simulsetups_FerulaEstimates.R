## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Model Selection Performance in Phylogenetic Comparative Methods under multivariate Ornstein–Uhlenbeck Models of Trait Evolution"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

load(paste0(ferula_estimation_directory,"BestModels_FerulaEstimation.RData"))

lmodelcomparison_setups<-sapply(1:8,function(i,lFerulaEstimPars){
    list(
	id=paste0("0",i),
	simulsetup=list(
	    evolmodel="ouch",
	    Atype= "Any",
	    Syytype= "Diagonal",
	    N=c(128,256,512,1024,2048),
            iter=list("size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100),
	    modelparams=lFerulaEstimPars[[i]],
	    Merror=NULL
	),
	recoverunder=list(
	    model_setups=list(
		list(evolmodel="bm"),
		list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c(NA,NA,"0","0","0"),c(NA,NA,"0","0","0"),c("0","0",NA,NA,"0"),c("0","0",NA,NA,"0"),c("0","0","0","0",NA)))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c(NA,NA,"0","0","0"),c(NA,NA,"0","0","0"),c("0","0",NA,NA,"0"),c("0","0",NA,NA,"0"),c("0","0","0","0",NA)))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c(NA,"0","0","0","0"),c("0",NA,"0","0","0"),c("0","0",NA,"0","0"),c("0","0","0",NA,"0"),c("0","0","0","0",NA)))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),
    		list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c(NA,"0","0","0","0"),c("0",NA,"0","0","0"),c("0","0",NA,"0","0"),c("0","0","0",NA,"0"),c("0","0","0","0",NA))))
	    ),
	    repeats=3,
	    kY=NULL,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    model_abbrev=c("bm_noME","m4SD_noME","m4SD_noMEdANULL","m4UT_noME","m4UT_noMEdANULL","m7SD_noME","m7SD_noMEdANULL","m7UT_noME","m7UT_noMEdANULL")
	)
    )},lFerulaEstimPars=lFerulaEstimPars,simplify=FALSE
)

lmodelcomparison_setups[[9]]<-
    list(
	id="09",
	simulsetup=list(
	    evolmodel="bm",
	    Atype= NULL,
	    Syytype= NULL,
	    N=c(128,256,512,1024,2048),
            iter=list("size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100),
	    modelparams=lFerulaEstimPars[[9]],
	    Merror=NULL
	),
	recoverunder=list(
	    model_setups=list(
		list(evolmodel="bm"),
		list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c(NA,NA,"0","0","0"),c(NA,NA,"0","0","0"),c("0","0",NA,NA,"0"),c("0","0",NA,NA,"0"),c("0","0","0","0",NA)))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c(NA,NA,"0","0","0"),c(NA,NA,"0","0","0"),c("0","0",NA,NA,"0"),c("0","0",NA,NA,"0"),c("0","0","0","0",NA)))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c(NA,"0","0","0","0"),c("0",NA,"0","0","0"),c("0","0",NA,"0","0"),c("0","0","0",NA,"0"),c("0","0","0","0",NA)))),
		list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),
    		list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c(NA,"0","0","0","0"),c("0",NA,"0","0","0"),c("0","0",NA,"0","0"),c("0","0","0",NA,"0"),c("0","0","0","0",NA))))
	    ),
	    repeats=3,
	    kY=NULL,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    model_abbrev=c("bm_noME","m4SD_noME","m4SD_noMEdANULL","m4UT_noME","m4UT_noMEdANULL","m7SD_noME","m7SD_noMEdANULL","m7UT_noME","m7UT_noMEdANULL")
	)
    )
lmodelcomparison_setups[[10]]<-lmodelcomparison_setups[[9]]
lmodelcomparison_setups[[10]]$simulsetup$modelparams<-lFerulaEstimPars[[10]]
lmodelcomparison_setups[[10]]$id<-"10"

v_df_invchisq<-rep(NA,5)
for (i in 1:5){## we have 5D data
    v_df_invchisq[i]<-1/mean(mFerula_Merror[,i],na.rm=TRUE)+2 ## mean of invchisq(df) is E[X]=1/(df-2)    
}
f_Merror<-function(N,v_df){
    Merror<-NULL
    ## inv-chisq(df) is 1/chisq(df)
    f_make_merror_matrix<-function(i,v_df){diag(1/(sapply(v_df,function(df_value){res<-NA;if(!is.na(df_value)){res<-rchisq(1,df_value)};res},simplify=TRUE)))}
    if (N==1){Merror<-f_make_merror_matrix(1,v_df)}else{Merror<-sapply(1:N,f_make_merror_matrix,v_df=v_df,simplify=FALSE)}
    Merror
}

for (i in seq(2,10,2)){
## even models are with measurement error
    lmodelcomparison_setups[[i]]$simulsetup$Merror<-f_Merror
    lmodelcomparison_setups[[i]]$simulsetup$Merror_parameters<-v_df_invchisq
    lmodelcomparison_setups[[i]]$recoverunder$model_abbrev<-c("bm_ME","m4SD_ME","m4SD_MEdANULL","m4UT_ME","m4UT_MEdANULL","m7SD_ME","m7SD_MEdANULL","m7UT_ME","m7UT_MEdANULL")
}

for (i in c(3,4,7,8)){
    lmodelcomparison_setups[[i]]$simulsetup$Syytype<-"UpperTri"
}


