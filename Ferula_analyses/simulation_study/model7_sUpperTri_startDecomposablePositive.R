## This file accompanies the manuscript: 
## Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Model Selection Performance in Phylogenetic Comparative Methods under multivariate Ornstein–Uhlenbeck Models of Trait Evolution"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

load("Ferula_Data.RData")
seed1 <- 10+541
set.seed(seed1)
OUOUmodel_start <- mvSLOUCH::ouchModel(ferulatreeape, feruladata_noME, regimes=NULL, Atype="DecomposablePositive", Syytype="UpperTri", diagA=NULL, maxiter=c(50,100))
seed2 <- 1991+541
set.seed(seed2)
OUOUmodel <- mvSLOUCH::ouchModel(ferulatreeape, feruladata_noME, regimes=NULL, Atype="Any", Syytype="UpperTri", diagA=NULL, maxiter=c(100,1000), start_point_for_optim=list(A=OUOUmodel_start$FinalFound$ParamsInModel$A, Syy=OUOUmodel_start$FinalFound$ParamsInModel$Syy), parameter_signs=list(signsA=rbind(c("+","0","0","0","0"),c("0","+","0","0","0"),c("0","0","+","0","0"),c("0","0","0","+","0"),c("0","0","0","0","+"))))
