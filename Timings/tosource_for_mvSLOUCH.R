## This file accompanies the manuscript: 
## Bartoszek, Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(mvSLOUCH) ## v 2.7.3
library(mvSLOUCHold) ## v 1.3.4, needs to be installed manually
library(TreeSim)
library(ouch)

f_loadsave_randomseed<-function(b_use_random_seed_from_manuscript,b_should_random_seed_be_saved,results_directory,filename_prefix,filename_RandomSeed,filename_suffix){
    rexp(1)
    crandomseed_filename<-paste0(results_directory,filename_prefix,filename_RandomSeed,filename_suffix,".RData")
    if(b_use_random_seed_from_manuscript){##setup the random number generator if replicating results
        b_randomseed_file_exists<-FALSE
        if(file.exists(crandomseed_filename)){		    
	    load(crandomseed_filename)
        Rseed<-.Random.seed
        rm(.Random.seed)
        RNGkind(kind = RNG_kind[1], normal.kind = RNG_kind[2], sample.kind = RNG_kind[3])
        RNGversion(RNG_version)     
	    assign('.Random.seed', Rseed, pos=.GlobalEnv)
        b_randomseed_file_exists<-TRUE
    }
    }
    if((!b_use_random_seed_from_manuscript)||(!b_randomseed_file_exists)){
	if (b_should_random_seed_be_saved){
    	    RNG_kind<-RNGkind()
    	    RNG_version<-getRversion()
    	    save(.Random.seed,RNG_kind,RNG_version,file=crandomseed_filename)
	}
    }	
}
