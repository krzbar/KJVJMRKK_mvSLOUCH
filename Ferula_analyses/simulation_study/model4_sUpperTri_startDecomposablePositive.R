load("Ferula_Data.RData")
seed1 <- 10+541 
set.seed(seed1)
OUOUmodel_start <- mvSLOUCH::ouchModel(ferulatreeape, feruladata_noME, regimes=NULL, Atype="DecomposablePositive", Syytype="UpperTri", diagA=NULL, maxiter=c(50,100))
seed2 <- 1991+541
set.seed(seed2)
OUOUmodel <- mvSLOUCH::ouchModel(ferulatreeape, feruladata_noME, regimes=NULL, Atype="Any", Syytype="UpperTri", diagA=NULL, maxiter=c(100,1000), start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+"))))


