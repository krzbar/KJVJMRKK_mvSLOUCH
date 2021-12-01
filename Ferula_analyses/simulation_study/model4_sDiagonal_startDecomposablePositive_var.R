load("Ferula_Data.RData")
seed1 <- 10+226
set.seed(seed1)
OUOUmodel_start <- mvSLOUCH::ouchModel(ferulatreeape, feruladata_ME, regimes=NULL, Atype="DecomposablePositive", Syytype="Diagonal", diagA=NULL, maxiter=c(25,50))
seed2 <- 1991+226
set.seed(seed2)
OUOUmodel <- mvSLOUCH::ouchModel(ferulatreeape, feruladata_ME, regimes=NULL, Atype="Any", Syytype="Diagonal", diagA=NULL, maxiter=c(100,1000), M.error=M.error, start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), parameter_signs=list(signsA=rbind(c("+",NA,"0","0","0"),c(NA,"+","0","0","0"),c("0","0","+",NA,"0"),c("0","0",NA,"+","0"),c("0","0","0","0","+"))))
