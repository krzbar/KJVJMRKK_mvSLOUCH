load("Carnivora_mvSLOUCH_objects_Full.RData") ## File with all estimation results by mvSLOUCH
BT$paramatric.bootstrap.estimation.replicates<-NULL
FinalOUf1$FinalFound$ParamSummary<-NULL
FinalOUf2$FinalFound$ParamSummary <- list(aic.c=FinalOUf2$FinalFound$ParamSummary$aic.c)
FinalOUs1$FinalFound$ParamSummary<-NULL
## nothing is removed from FinalOUs2 as substantial parts of ParamSummary are required
tmpOU1testedmodels<-OU1$testedModels
OU1<-list(model.setups=OU1$model.setups,BestModel=OU1$BestModel,testedModels=as.list(rep(NA,length(OU1$testedModels))))
OU1$testedModels[[11]]<-tmpOU1testedmodels[[11]]
OU1$testedModels[[21]]<-tmpOU1testedmodels[[21]]
rm("tmpOU1testedmodels")
OU1BM$FinalFound$ParamSummary<-list(phyl.halflife=OU1BM$FinalFound$ParamSummary$phyl.halflife,optimal.regression=OU1BM$FinalFound$ParamSummary$optimal.regression,evolutionary.regression=OU1BM$FinalFound$ParamSummary$evolutionary.regression,RSS=OU1BM$FinalFound$ParamSummary$RSS,dof=OU1BM$FinalFound$ParamSummary$dof)
OU1OU$FinalFound$ParamSummary<-list(phyl.halflife=OU1OU$FinalFound$ParamSummary$phyl.halflife,RSS=OU1OU$FinalFound$ParamSummary$RSS,evolutionary.regression=OU1OU$FinalFound$ParamSummary$evolutionary.regression,dof=OU1OU$FinalFound$ParamSummary$dof)
OUBMestim ## all is needed
OUBMestim.mod<-list(MaxLikFound=OUBMestim.mod$MaxLikFound)
OUOUestim$FinalFound$ParamSummary<-list(aic.c=OUOUestim$FinalFound$ParamSummary$aic.c,phyl.halflife=OUOUestim$FinalFound$ParamSummary$phyl.halflife,RSS=OUOUestim$FinalFound$ParamSummary$RSS,evolutionary.regression=OUOUestim$FinalFound$ParamSummary$evolutionary.regression,dof=OUOUestim$FinalFound$ParamSummary$dof)
OUOUreStart$FinalFound$ParamSummary<-NULL
OUOUstart$FinalFound$ParamSummary<-NULL
OUc<-list(BestModel=list(model=OUc$BestModel$model,BestModel=OUc$BestModel$BestModel))
OUf<-list(BestModel=list(model=OUf$BestModel$model,BestModel=OUf$BestModel$BestModel))
OUr<-list(BestModel=list(model=OUr$BestModel$model,BestModel=OUr$BestModel$BestModel))
OUs<-list(BestModel=list(model=OUs$BestModel$model,BestModel=OUs$BestModel$BestModel))
OptOUf1$FinalFound<-list(ParamsInModel=OptOUf1$FinalFound$ParamsInModel,LogLik=OptOUf1$FinalFound$LogLik)
OptOUf2$FinalFound<-list(LogLik=OptOUf2$FinalFound$LogLik)
OptOUs1$FinalFound<-list(ParamsInModel=OptOUs1$FinalFound$ParamsInModel,LogLik=OptOUs1$FinalFound$LogLik)
OptOUs2$FinalFound<-list(LogLik=OptOUs2$FinalFound$LogLik)
save.image(file = "Carnivora_mvSLOUCH_objects.RData") ## RData file that is included with mvSLOUCH