load("Ferula_Data.RData")
seed1 <- 10+90
set.seed(seed1)
OUOUmodel_start <- mvSLOUCH::ouchModel(ferulatreeape, feruladata_noME, regimes=NULL, Atype="DecomposablePositive", Syytype="Diagonal", diagA=NULL, maxiter=c(50,100))
seed2 <- 1991+90
set.seed(seed2)
OUOUmodel <- mvSLOUCH::ouchModel(ferulatreeape, feruladata_noME, regimes=NULL, Atype="Any", Syytype="Diagonal", diagA=NULL, maxiter=c(100,1000), start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+"))))
