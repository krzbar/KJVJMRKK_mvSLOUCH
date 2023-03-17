## This file accompanies the manuscript: Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje "Analytical advances alleviate model misspecification in non--Brownian multivariate comparative methods"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

lmodelcomparison_setups<-sapply(1:length(v_p),function(i,vN,viter,v_p){
    p<-v_p[i]
    id_char<-as.character(i)
    if (i<10){id_char<-paste0("0",i)}
    list(
	id=id_char,
	simulsetup=list(
	    evolmodel="bm_indep",
	    Atype= NULL,
	    Syytype= NULL,
	    ## these are now defined in the usercontrolledvariables.R file
#	    N=c(32,64,128,256,512,1024,2048),
#            iter=list("size_32"=1:100,"size_64"=1:100,"size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100),
	    N=vN,
            iter=viter,
	    modelparams=list(vX0=matrix(rep(0,p),ncol=1),Sxx=diag(1,p,p)),
	    Merror=NULL
	),
	recoverunder=list(
	    model_setups=list(
		list(evolmodel="bm"),
		list(evolmodel="ouch","Atype"="SingleValueDiagonal","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="SingleValueDiagonal","Syytype"="Diagonal","diagA"=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL),
		list(evolmodel="mvslouch","Atype"="SingleValueDiagonal","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="mvslouch","Atype"="SingleValueDiagonal","Syytype"="Diagonal","diagA"=NULL),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL)
	    ),
	    non_mvSLOUCH_model_setups=list(
	    ## additional models that are not present in mvSLOUCH
	    	list(evolmodel="bm_indep")
	    ),
	    repeats=3,
	    kY=2,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    ##model_abbrev=c("bm","bm_indep","OUOU1DP","OUOU1D","OUOU4DP","OUOU4D","OUBM1DPML","OUBM1DML","OUBM4DPML","OUBM4DML")
	    ##model_abbrev=c("bm","bm_indep","OUOU1DP","OUOU1D","OUOU4DP","OUOU4D","OUBM1DPML","OUBM1DML","OUBM4DPML","OUBM4DML")
	    model_abbrev=c("bm","OUOU1DP","OUOU1D","OUOU4DP","OUOU4D","OUBM1DPML","OUBM1DML","OUBM4DPML","OUBM4DML","bm_indep")
	)
    )},vN=vN,viter=viter,v_p=v_p,simplify=FALSE
)

