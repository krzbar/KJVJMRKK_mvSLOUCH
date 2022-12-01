## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

load(file=paste0("./",c_resultsdirectory,"/initialdiagUTASR/regimes_Ellenberg_tree.RData"))


## MODEL BM  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
result<-multiResultClass()
if (b_use_random_seed_from_manuscript){
        load(paste0("./",c_rseedoutfiledirectoryprefix,"/modelBM/",c_rseedoutfiledirectoryprefix,"_modelBM.RData"))
        RNGkind(kind = RNG_kind_modelBM[1], normal.kind = RNG_kind_modelBM[2], sample.kind = RNG_kind_modelBM[3]) 
        RNGversion(min(as.character(getRversion(),RNG_version_modelBM)))
        rm(.Random.seed)
        assign('.Random.seed', random_seed_modelBM, pos=.GlobalEnv)
}else{
        RNG_kind_modelBM<-RNGkind()
        RNG_version_modelBM<-getRversion()
        random_seed_modelBM<-.Random.seed
}
if (b_save_random_seeds){save(random_seed_modelBM,RNG_kind_modelBM,RNG_version_modelBM,file=paste0("./",c_rseedoutfiledirectoryprefix,"/modelBM/",c_rseedoutfiledirectoryprefix,"_modelBM.RData"))}

start_time <- Sys.time()

result$mvSLOUCH_output <- mvSLOUCH::BrownianMotionModel(mvStree, mvData2)
  
end_time <- Sys.time(); end_time - start_time 


result$runtime <- end_time - start_time 
result$random_seed <-random_seed_modelBM
res_BM<-result$mvSLOUCH_output;runtime<-result$runtime;save(res_BM, runtime,file=paste0("./",c_indivresultsdirectory ,"/modelBM/OUOU.modelBM.RData"))
BM.model<-list()
BM.model[[1]]<-result
saveRDS(BM.model, paste0("./",c_resultsdirectory,"/modelBM/BMmodel.rds"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
