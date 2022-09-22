## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Model Selection Performance in Phylogenetic Comparative Methods under multivariate Ornstein–Uhlenbeck Models of Trait Evolution"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


N<-length(ferulatreeape$tip.label)
lmodelcomparison_setups<-sapply(1:8,function(i,ferulatreeape,lFerulaEstimPars){
    list(
	id=paste0("0",i),
	simulsetup=list(
	    evolmodel="ouch",
	    Atype= "Any",
	    Syytype= "Diagonal",
	    N=c(N),
	    iter=list(size_N=1:100),
	    modelparams=lFerulaEstimPars[[i]],
	    ape_phyltree=ferulatreeape,
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
		##list(evolmodel="bm"),
		##list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),
		##list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+")))),
		##list(evolmodel="ouch","Atype"="Any","Syytype"="Diagonal","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+")))),
		##list(evolmodel="ouch","Atype"="Any","Syytype"="UpperTri","diagA"=NULL,parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+"))))
	    ),
	    repeats=3,
	    ##kY=5,
	    kY=NULL,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    model_abbrev=c("bm_noME","m4SD_noME","m4SD_noMEdANULL","m4UT_noME","m4UT_noMEdANULL","m7SD_noME","m7SD_noMEdANULL","m7UT_noME","m7UT_noMEdANULL")
	    ##c("bm_noME","m4SD_noME","m4UT_noME","m7SD_noME","m7UT_noME")
	)
    )},ferulatreeape=ferulatreeape,lFerulaEstimPars=lFerulaEstimPars,simplify=FALSE
)

lmodelcomparison_setups[[9]]<-
    list(
        id="09",
        simulsetup=list(
            evolmodel="bm",
            Atype= NULL,
            Syytype= NULL,
            N=c(N),
            iter=list(size_N=1:100),
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
            model_abbrev= c("bm_noME","m4SD_noME","m4SD_noMEdANULL","m4UT_noME","m4UT_noMEdANULL","m7SD_noME","m7SD_noMEdANULL","m7UT_noME","m7UT_noMEdANULL")
        )
    )
lmodelcomparison_setups[[10]]<-lmodelcomparison_setups[[9]]
lmodelcomparison_setups[[10]]$simulsetup$modelparams<-lFerulaEstimPars[[10]]
lmodelcomparison_setups[[10]]$id<-"10"

for (i in seq(2,10,2)){
## even models are with measurement error
    lmodelcomparison_setups[[i]]$simulsetup$Merror<-M.error
    lmodelcomparison_setups[[i]]$recoverunder$model_abbrev<-c("bm_ME","m4SD_ME","m4SD_MEdANULL","m4UT_ME","m4UT_MEdANULL","m7SD_ME","m7SD_MEdANULL","m7UT_ME","m7UT_MEdANULL")
    ##c("bm_ME","m4SD_ME","m4UT_ME","m7SD_ME","m7UT_ME")
}

for (i in c(3,4,7,8)){
    lmodelcomparison_setups[[i]]$simulsetup$Syytype<-"UpperTri"
}

lmodelcomparison_setups<-sapply(lmodelcomparison_setups,function(lsetup,N){
    names(lsetup$simulsetup$iter)<-paste0("size_",N)
    lsetup
},N=N,simplify=FALSE)

