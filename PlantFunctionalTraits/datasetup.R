## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

RNG_kind<-RNGkind()
RNG_version<-getRversion()


## =============== Read files ============================================================================
#Tree (scaled to unit length)
ScaledTree<-read.tree(paste0("./",c_datadirectory,"/",c_phylogeny_file))

# plant data (including data for computing regimes)
dat<-read.csv(paste0("./",c_datadirectory,"/",c_trait_measurements_file), row.names = 1, header=TRUE)
## =============== End read files ========================================================================


## ====================== Trait data preparation ===============================================================
### Left untouched in case somebody (for whatever reason) wants to have access to the data on the original scale.            

# Before running analyses, we will log transform the morphological variables so that they are less 
# susceptible to scaling effects:
mvData <- data.matrix(dat[ ,c("leaf.area", "plant.height", "seed.mass", "leaf.mass")])
mvData <- log(mvData)
mvData2 <- mvData[ , c("plant.height", "seed.mass", "leaf.area", "leaf.mass")]

## ====================== End of trait data preparation ========================================================



## ======= Create required directory structure ===========================================================
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/model1/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/model2/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/model3/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/model4/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/model5/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/model6/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/bootstrap/"), showWarnings = FALSE)
dir.create(paste0("./",c_rseedoutfiledirectoryprefix,"/modelBM/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/initialdiagUTASR/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/model1/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/model2/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/model3/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/model4/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/model5/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/model6/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/bootstrap/"), showWarnings = FALSE)
dir.create(paste0("./",c_resultsdirectory,"/modelBM/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/initialdiagUTASR/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/model1/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/model2/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/model3/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/model4/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/model5/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/model6/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/bootstrap/"), showWarnings = FALSE)
dir.create(paste0("./",c_indivresultsdirectory,"/modelBM/"), showWarnings = FALSE)
## ======= End create required directory structure =======================================================

## ===== Initialize random seeds =========================================================================
rexp(1) ## This is important as it will initialize the random number generator
## read in the random seed information
## the below file needs to have the variables RNG_kind, RNG_version and random seed of ASR (random_seed_ASR) and the initial diagonal run
if ((b_use_random_seed_from_manuscript) || (b_use_random_seed_from_manuscript_forASR)){
    load(paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASR_initial.RData"))
    RNGkind(kind = RNG_kind[1], normal.kind = RNG_kind[2], sample.kind = RNG_kind[3]) 
    RNGversion(min(as.character(getRversion(),RNG_version)))
    assign('.Random.seed', random_seed_initial, pos=.GlobalEnv)
}else{
    RNG_kind<-RNGkind()
    RNG_version<-getRversion()
    random_seed_initial<-.Random.seed
    save(RNG_kind,random_seed_initial,RNG_version,file=paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASR_initial.RData"))
}
## ===== End initialize random seeds =====================================================================
