## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

# Set number of cores
num_cores <- 48

### Use multiple cores: TRUE/FALSE
run.parallel<-TRUE 

### Rerun original code: TRUE/FALSE. If FALSE, a new set of seeds is used (and saved if b_save_random_seeds is TRUE).
b_use_random_seed_from_manuscript <-FALSE
b_save_random_seeds<-TRUE ## TRUE/FALSE should the random seeds be saved
b_use_random_seed_from_manuscript_forASR <-FALSE  ## if we want to force only the regime calculations to be the same

c_rseedoutfiledirectoryprefix<-"RandomSeeds"
c_resultsdirectory<- "Results"
c_indivresultsdirectory<- "IndivResults"

c_datadirectory<- "Data"
c_trait_measurements_file<-"plant_data.csv"
c_phylogeny_file<-"plant_tree.txt"

b_doMLresolution_amb <- FALSE ##TRUE/FALSE if TRUE for ambiguous parsimony based regimes the regime value is taken from the simulated posterior mode, if FALSE manual resolution is done
num_trees_regimes_simul <- 5000 ## number of simulation runs to get probability of given regime for each branch for phytools::make.simmap()

num_model_repeats<- 500 ## perhaps make less for code testing
num_bootstrap_repeats<- 500 ## perhaps make less for code testing
#num_model_repeats<- 3 ## perhaps make less for code testing
#num_bootstrap_repeats<- 3 ## perhaps make less for code testing

v_statistics_to_bootstrap<-c("corr.matrix", "trait.regression", "phyl.halflife") ## statistics for which to create the bootstrap CIs
bootci_lvl<-0.95 ## what is the size of the bootstrap CIs, here 95%

b_make_optima_plot <-TRUE ## should the plot with the optima values be made TRUE/FALSE
b_make_phyl_regimes_plot <- TRUE ## should the phylogeny with the regimes be plotted TRUE/FALSE

