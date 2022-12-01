## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

##==================== Regime creation ===================================================================
if (b_use_random_seed_from_manuscript||b_use_random_seed_from_manuscript_forASR){.Random.seed<-random_seed_ASR}else{random_seed_ASR<-.Random.seed}

## Jesualdo please edit below the code to create regimes without ambiguities 
## perhaps one could compare Fitch regimes with other ones for a short paragraph in the Supplementary Material if biologically or from a user's perspective interesting?

#Storing the Nutrient habitat categories in a new object where the order is specified according to tree tip names:
regimes_dat <- dat$Nutrients[order(match(row.names(dat), ScaledTree$tip.label))]
regimes_dat<-paste0("Elvl_",regimes_dat) ## mvSLOUCH::fitch.mvsl() can behave strange if the regimes are just numbers

# With this new object we can run the parsimony reconstruction:
regimesFitch_nodeltran <- mvSLOUCH::fitch.mvsl(ScaledTree, regimes_dat)

# A lot of ambigous nodes so we did a delayed transformations reconstruction (which "only" gives 19 ambigous nodes):
regimesFitch <- mvSLOUCH::fitch.mvsl(ScaledTree, regimes_dat, deltran=TRUE)

## save the regimes vector to be used
## Here we have those regimes that will be used, please choose the one that Jesualdo creates
regimes_EllenbergFitch<- regimesFitch$branch_regimes 
root_regime<- as.character(regimesFitch$root_regime)

regimes_FitchML<-NULL
mtreeER<-NULL
mtreeSYM<-NULL
mtreeARD<-NULL
mASRmodel<-NULL
c_ASR_model<-NULL
mtrees<-NULL
XX<-NULL
regimes_ML<-NULL
aicER<-NA
aiccER<-NA
aicSYM<-NA
aiccSYM<-NA
aicARD<-NA
aiccARD<-NA

if (b_doMLresolution_amb){
## =========================== ML regimes creation =====================================================
    # Stochastic character mapping (MP-informed root)
    regimes_dat<-as.character(regimes_dat)
    names(regimes_dat)<-ScaledTree$tip.label
    print(geiger::name.check(ScaledTree,regimes_dat))

aic<-function(logL,k) 2*k-2*logL
aicc<-function(aic,k,N) aic+((2*k*(k+1))/(N-k-1))

    # 1. Equal-rates
    if ((b_use_random_seed_from_manuscript)||(b_use_random_seed_from_manuscript_forASR)){ 
        	load(paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASRER.RData"))
	    RNGkind(kind = RNG_kind_ASRER[1], normal.kind = RNG_kind_ASRER[2], sample.kind = RNG_kind_ASRER[3]) 
    	    RNGversion(min(as.character(getRversion(),RNG_version_ASRER)))
        	rm(.Random.seed)
	    assign('.Random.seed', random_seed_ASRER, pos=.GlobalEnv)
        }else{
	    RNG_kind_ASRER<-RNGkind()
            RNG_version_ASRER<-getRversion()
    	    random_seed_ASRER<-.Random.seed
	}
    if (b_save_random_seeds){save(random_seed_ASRER,RNG_kind_ASRER,RNG_version_ASRER,file=paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASRER.RData"))}
    mtreeER <- phytools::make.simmap(ScaledTree, regimes_dat, Q="empirical", model = "ER", pi = "equal")
    aicER<-aic(mtreeER$logL,1)
    aiccER<-aicc(aicER,1,length(ScaledTree$tip.label))

    # 2. Symmetrical
    if ((b_use_random_seed_from_manuscript)||(b_use_random_seed_from_manuscript_forASR)){ 
	    load(paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASRSYM.RData"))
            RNGkind(kind = RNG_kind_ASRSYM[1], normal.kind = RNG_kind_ASRSYM[2], sample.kind = RNG_kind_ASRSYM[3]) 
            RNGversion(min(as.character(getRversion(),RNG_version_ASRSYM)))
		rm(.Random.seed)
	    assign('.Random.seed', random_seed_ASRSYM, pos=.GlobalEnv)
        }else{
	    RNG_kind_ASRSYM<-RNGkind()
            RNG_version_ASRSYM<-getRversion()
    	    random_seed_ASRSYM<-.Random.seed
	}
    if (b_save_random_seeds){save(random_seed_ASRSYM,RNG_kind_ASRSYM,RNG_version_ASRSYM,file=paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASRSYM.RData"))}

    mtreeSYM <- phytools::make.simmap(ScaledTree, regimes_dat, Q="empirical", model = "SYM", pi = "equal")
    aicSYM<-aic(mtreeSYM$logL,sum(upper.tri(mtreeSYM$Q)))
    aiccSYM<-aicc(aicSYM,sum(upper.tri(mtreeSYM$Q)),length(ScaledTree$tip.label))

    # 3. All-rates-different
    if ((b_use_random_seed_from_manuscript)||(b_use_random_seed_from_manuscript_forASR)){ 
	    load(paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASRARD.RData"))
            RNGkind(kind = RNG_kind_ASRARD[1], normal.kind = RNG_kind_ASRARD[2], sample.kind = RNG_kind_ASRARD[3]) 
	    RNGversion(min(as.character(getRversion(),RNG_version_ASRARD)))
            rm(.Random.seed)
	    assign('.Random.seed', random_seed_ASRARD, pos=.GlobalEnv)
        }else{
	    RNG_kind_ASRARD<-RNGkind()
            RNG_version_ASRARD<-getRversion()
    	    random_seed_ASRARD<-.Random.seed
	}
    if (b_save_random_seeds){save(random_seed_ASRARD,RNG_kind_ASRARD,RNG_version_ASRARD,file=paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASRARD.RData"))}

    mtreeARD <- phytools::make.simmap(ScaledTree, regimes_dat, Q="empirical", model = "ARD", pi = "equal")
    aicARD<-aic(mtreeARD$logL,sum(upper.tri(mtreeARD$Q))*2)
    aiccARD<-aicc(aicARD,sum(upper.tri(mtreeARD$Q))*2,length(ScaledTree$tip.label))

    index_min_aicc<-which.min(c(aiccER,aiccSYM,aiccARD))
    mASRmodel<-switch(index_min_aicc,mtreeER,mtreeSYM,mtreeARD)

    c_ASR_model<-switch(index_min_aicc,"ER","SYM","ARD")

    if ((b_use_random_seed_from_manuscript)||(b_use_random_seed_from_manuscript_forASR)){ 
    	    load(paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASRsimbest.RData"))
	    RNGkind(kind = RNG_kind_ASRsimbest[1], normal.kind = RNG_kind_ASRsimbest[2], sample.kind = RNG_kind_ASRsimbest[3]) 
            RNGversion(min(as.character(getRversion(),RNG_version_ASRsimbest)))
    	    rm(.Random.seed)
	    assign('.Random.seed', random_seed_ASRsimbest, pos=.GlobalEnv)
        }else{
	    RNG_kind_ASRsimbest<-RNGkind()
            RNG_version_ASRsimbest<-getRversion()
    	    random_seed_ASRsimbest<-.Random.seed
	}
    if (b_save_random_seeds){save(random_seed_ASRsimbest,RNG_kind_ASRsimbest,RNG_version_ASRsimbest,file=paste0("./",c_rseedoutfiledirectoryprefix,"/initialdiagUTASR/",c_rseedoutfiledirectoryprefix,"_ASRsimbest.RData"))}

    mtrees <- phytools::make.simmap(ScaledTree, regimes_dat, Q="empirical", model = c_ASR_model, pi = "estimated", nsim = num_trees_regimes_simul)
    XX <- phytools::describe.simmap(mtrees)
    m_reg_postprobs<-XX$ace
    regimes_ML<-rep(NA,length(regimes_dat))

    for (i in 1:length(regimes_dat)){
	tip_index<-which(ScaledTree$tip.label==names(regimes_dat)[i])
        edge_index<-which(ScaledTree$edge[,2]==tip_index)
	regimes_ML[edge_index]<-regimes_dat[i]
    }
    for (i in 1:nrow(m_reg_postprobs)){
        c_ML_edge_regime<-colnames(m_reg_postprobs)[which.max(m_reg_postprobs[i,])]
	node_index<-as.numeric(row.names(m_reg_postprobs)[i])
        if (node_index==length(ScaledTree$tip.label)+1){## we are in the root regime
    	    root_regime<-c_ML_edge_regime
	}else{
            edge_index<-which(ScaledTree$edge[,2]==node_index)
	    regimes_ML[edge_index]<-c_ML_edge_regime
	}
    }
    regimes_FitchML<-regimes_EllenbergFitch
    v_which_amb<-which(regimes_EllenbergFitch=="ambiguous")
    if (length(v_which_amb)>0){
	regimes_FitchML[v_which_amb]<-regimes_ML[v_which_amb]
    }

    print(paste0("Number of positions at which Fitch and ML regimes differ: ",length(which(regimes_FitchML!=regimes_ML))))
    regimes_Ellenberg<-regimes_FitchML
## =========================== end of ML regimes creation ==============================================
}else{
    if (root_regime=="ambiguous"){
        root_regime<- names(table(regimes_EllenbergFitch))[which.max(table(regimes_EllenbergFitch))]
    }
    regimes_EllenbergFitch[791]<-"Elvl_5" ## clade with 1 ambiguous
    regimes_EllenbergFitch[c(247,256,258,259,260)]<-"Elvl_6" ## clade with 5 ambiguous
    ## clade with 13 ambiguous
    regimes_EllenbergFitch[940]<-"Elvl_3" 
    regimes_EllenbergFitch[937]<-"Elvl_4" 
    regimes_EllenbergFitch[c(926,927,928,930,931,932,945)]<-"Elvl_5" 
    regimes_EllenbergFitch[c(924)]<-"Elvl_6" 
    regimes_EllenbergFitch[c(947,948,949)]<-"Elvl_7" 
    ## ==========================================
    regimes_Ellenberg<-regimes_EllenbergFitch
}

regimes_Ellenberg_forplot<-regimes_Ellenberg
sorted_regs<-sort(unique(regimes_Ellenberg_forplot))
regs_cols<-setNames(c(
"blue",
"cadetblue4",
"cadetblue1",
"lightcyan",
"lavenderblush",
"bisque",
"chocolate1",
"red",
"red4"
), sorted_regs)

if (b_make_phyl_regimes_plot){
    regimes_Ellenberg_forplot_col<-rep(NA,length(regimes_Ellenberg_forplot))
    for(i in 1:length(regimes_Ellenberg_forplot_col)){
	regimes_Ellenberg_forplot_col[i]<-regs_cols[which(sorted_regs==regimes_Ellenberg_forplot[i])]
    }

    pdf(paste0("./",c_resultsdirectory,"/FigTreeRegimes.pdf"))
    plot(ScaledTree,cex = 0.3,edge.color = regimes_Ellenberg_forplot_col,type="fan",show.tip.label=FALSE)
    #legend("bottom",legend=as.character(1:length(regs_cols)),col=regs_cols[1:length(regs_cols)],bty="n",lty=1,lwd=3,ncol=5)
    #legend(-1,-0.95,legend=as.character(1:length(regs_cols)),col=regs_cols[1:length(regs_cols)],bty="n",lty=1,lwd=3,ncol=5)
    legend(-1,1.125,legend=as.character(1:5),col=regs_cols[1:5],bty="n",lty=1,lwd=3,ncol=5)
    legend(-1,-0.95,legend=as.character(6:9),col=regs_cols[6:9],bty="n",lty=1,lwd=3,ncol=5)
    dev.off()
}
##==================== End of regime creation ============================================================


# mvSLOUCH works faster when the phylogeny object includes information on the paths and distances 
# for each node. This information can be obtained with the function mvSLOUCH::phyltree_paths():
mvStree <- mvSLOUCH::phyltree_paths(ScaledTree) 
save(regimesFitch,regimesFitch_nodeltran,regimes_Ellenberg,ScaledTree,mvStree,regimes_FitchML,regimes_EllenbergFitch,mtreeER,mtreeSYM,mtreeARD,mASRmodel,c_ASR_model,mtrees,XX,regimes_ML,aicER,aiccER,aicSYM,aiccSYM,aicARD,aiccARD,file=paste0("./",c_resultsdirectory,"/initialdiagUTASR/regimes_Ellenberg_tree.RData"))

## ====================== End of data preparation ========================================================

