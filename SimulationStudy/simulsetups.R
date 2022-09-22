## This file accompanies the manuscripts: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Model Selection Performance in Phylogenetic Comparative Methods under multivariate Ornstein–Uhlenbeck Models of Trait Evolution"
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Rotation Invariance in non–Brownian motion Phylogenetic Comparative Methods: A comment on Adams and Collyer (2018)"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

lmodelcomparison_setups<-list(
    list(
	id="01",
	simulsetup=list(
	    evolmodel="bm",
	    Atype=NULL,
	    Syytype=NULL,
	    N=c(32,64,128,256,512,2048),
	    iter=list("size_32"=1:100,"size_64"=1:100,"size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100),
	    modelparams=list(vX0=matrix(0,nrow=4,ncol=1), Sxx=rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)))
	),
	recoverunder=list(
	    model_setups=list(
		list(evolmodel="bm"),
		list(evolmodel="ouch","Atype"="SingleValueDiagonal","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri",diagA=NULL),
		list(evolmodel="ouch","Atype"="Invertible","Syytype"="UpperTri","diagA"="Positive","OUOU5"),
		list(evolmodel="mvslouch","Atype"="SingleValueDiagonal","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"=NULL,estimateBmethod="ML")
	    ),
	    repeats=3,
	    kY=2,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    model_abbrev=c("bm","OUOU1","OUOU4P","OUOU4","OUOU5","OUBM1ML","OUBM4PML","OUBM4ML")
	)
    ),
    list(
	id="02",
	simulsetup=list(
	    evolmodel="ouch",
	    Atype="Diagonal",
	    Syytype="Diagonal",
	    N=c(32,64,128,256,512,1024,2048),
	    iter=list("size_32"=1:100,"size_64"=1:100,"size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100),
	    modelparams=list(vY0=matrix(c(0,0,0,0),nrow=4,ncol=1),A=rbind(c(1,0,0,0),c(0,2,0,0),c(0,0,3,0),c(0,0,0,4)),mPsi=matrix(c(0,0,0,0),nrow=4,ncol=1), Syy=rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)))
	),
	recoverunder=list(
	    model_setups=list(
		list(evolmodel="bm"),
		list(evolmodel="ouch","Atype"="Diagonal","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="Diagonal","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri",diagA=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="mvslouch","Atype"="Diagonal","Syytype"="Diagonal","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="Diagonal","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL,estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"=NULL,estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML")
	    ),
	    repeats=3,
	    kY=2,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    model_abbrev=c("bm","OUOU2","OUOU2U","OUOU3","OUOU3P","OUOU4","OUOU4P","OUBM2ML","OUBM2UML","OUBM3ML","OUBM3PML","OUBM4ML","OUBM4PML")
	)
    ),
    list(
	id="03",
	simulsetup=list(
	    evolmodel="ouch",
	    Atype="DecomposablePositive",
	    Syytype="Diagonal",
	    N=c(32,64,128,256,512,1024,2048),
	    iter=list("size_32"=1:100,"size_64"=1:100,"size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100),
	    modelparams=list(vY0=matrix(c(0,0,0,0),nrow=4,ncol=1),A=rbind(c(1,0.2,0.3,0.4),c(0.4,2,0.2,0.1),c(0.1,0.3,3,0.4),c(0.1,0.1,0.1,4)),mPsi=matrix(c(0,0,0,0),nrow=4,ncol=1), Syy=rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)))
	),
	recoverunder=list(
	    model_setups=list(
		list(evolmodel="bm"),
		list(evolmodel="ouch","Atype"="Diagonal","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="Diagonal","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri",diagA=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="mvslouch","Atype"="Diagonal","Syytype"="Diagonal","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="Diagonal","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL,estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"=NULL,estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML")
	    ),
	    repeats=3,
	    kY=2,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    model_abbrev=c("bm","OUOU2","OUOU2U","OUOU3","OUOU3P","OUOU4","OUOU4P","OUBM2ML","OUBM2UML","OUBM3ML","OUBM3PML","OUBM4ML","OUBM4PML")
	)
    ),
    list(
	id="04",
	simulsetup=list(
	    evolmodel="mvslouch",
	    Atype="Diagonal",
	    Syytype="Diagonal",
	    N=c(32,64,128,256,512,1024,2048),	    
	    iter=list("size_32"=1:100,"size_64"=1:100,"size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100),
	    modelparams=list(vY0=matrix(c(0,0),nrow=2,ncol=1),A=rbind(c(1,0),c(0,2)),B=rbind(c(5,-5),c(-5,5)),mPsi=matrix(c(0,0),nrow=2,ncol=1),Syy=rbind(c(1,0),c(0,1)),vX0=matrix(c(0,0),nrow=2,ncol=1),Sxx=rbind(c(2,0),c(0,2)),Syx=matrix(0,ncol=2,nrow=2),Sxy=matrix(0,ncol=2,nrow=2))
	),
	recoverunder=list(
	    model_setups=list(
		list(evolmodel="bm"),
		list(evolmodel="ouch","Atype"="Diagonal","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="Diagonal","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri",diagA=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="mvslouch","Atype"="Diagonal","Syytype"="Diagonal","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="Diagonal","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL,estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"=NULL,estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML")
	    ),
	    repeats=3,
	    kY=2,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    model_abbrev=c("bm","OUOU2","OUOU2U","OUOU3","OUOU3P","OUOU4","OUOU4P","OUBM2ML","OUBM2UML","OUBM3ML","OUBM3PML","OUBM4ML","OUBM4PML")
	)
    ),
    list(
	id="05",
	simulsetup=list(
	    evolmodel="mvslouch",
	    Atype="DecomposablePositive",
	    Syytype="Diagonal",
	    N=c(32,64,128,256,512,1024,2048),
	    iter=list("size_32"=1:100,"size_64"=1:100,"size_128"=1:100,"size_256"=1:100,"size_512"=1:100,"size_1024"=1:100,"size_2048"=1:100),
	    modelparams=list(vY0=matrix(c(0,0),nrow=2,ncol=1),A=rbind(c(1,0.1),c(0.3,2)),B=rbind(c(5,-5),c(-5,5)),mPsi=matrix(c(0,0),nrow=2,ncol=1),Syy=rbind(c(1,0),c(0,1)),vX0=matrix(c(0,0),nrow=2,ncol=1),Sxx=rbind(c(2,0),c(0,2)),Syx=matrix(0,ncol=2,nrow=2),Sxy=matrix(0,ncol=2,nrow=2))
	),
	recoverunder=list(
	    model_setups=list(
		list(evolmodel="bm"),
		list(evolmodel="ouch","Atype"="Diagonal","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="Diagonal","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive"),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri",diagA=NULL),
		list(evolmodel="ouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive"),
		list(evolmodel="mvslouch","Atype"="Diagonal","Syytype"="Diagonal","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="Diagonal","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"=NULL,estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="Diagonal","diagA"="Positive",estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"=NULL,estimateBmethod="ML"),
		list(evolmodel="mvslouch","Atype"="DecomposablePositive","Syytype"="UpperTri","diagA"="Positive",estimateBmethod="ML")
	    ),
	    repeats=3,
	    kY=2,
	    ## abbreviations are for easiness of downstream table construction summarizing the model selection
	    model_abbrev=c("bm","OUOU2","OUOU2U","OUOU3","OUOU3P","OUOU4","OUOU4P","OUBM2ML","OUBM2UML","OUBM3ML","OUBM3PML","OUBM4ML","OUBM4PML")
	)
    )
)
