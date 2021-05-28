These are the R scripts and numerical results accompanying Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies".

The R setup for the manuscript was as follows: R version 3.6.1 (2019-09-12) Platform: x86_64-pc-linux-gnu (64-bit) Running under: openSUSE Leap 42.3

The exact output can depend on the random seed. However, in the script we have the option of rerunning the analyses as it was in the manuscript, i.e.
the random seeds that were used to generate the results are saved, included and can be read in.

The code is divided into several directories with scripts, random seeds and result files.

1) LikelihoodTesting

    Directory contains the script test_rotation_invariance_mvSLOUCH.R that demonstrates that mvSLOUCH's likelihood calculations are rotation invariant.
        
2) Carnivorans

    Directory contains files connected to the Carnivrons' vignette in mvSLOUCH.
    
    2.1) Carnivora_mvSLOUCH_objects_Full.RData
        Full output of  running the R code in the vignette. With mvSLOUCH is a very bare-minimum subset of this file that allows for the creation of the            vignette.
        
    2.2) Carnivora_mvSLOUCH_objects.RData
        Reduced objects from Carnivora_mvSLOUCH_objects_Full.RData that are included with mvSLOUCH's vignette.
        
    2.3) Carnivora_mvSLOUCH_objects_remove_script.R
        R script to reduce Carnivora_mvSLOUCH_objects_Full.RData to Carnivora_mvSLOUCH_objects.RData .
        
    2.4) mvSLOUCH_Carnivorans.Rmd
        The vignette itself.
        
    2.5) refs_mvSLOUCH.bib 
        Bib file for the vignette.
        
    2.6) ScaledTree.png, ScaledTree2.png, ScaledTree3.png, ScaledTree4.png
        Plots of phylogenetic trees for vignette.

3) SimulationStudy

    Directory contains all the output of the simulation study presented in the manuscript and scripts that allow for replication (the random number generator seeds are also provided) or running ones own simulation study, and scripts to generate graphs, and model comparison summary.
    
4) Ungulates

    Directory contains files connected to the "Feeding styles and oral morphology in ungulates" analyses performed for the manuscript.
    
    4.1) Data.csv
        The phenotypic data includes three continuous variables and one categorical variable. Continuous variables (MZW: muzzle width; HM3: unworn lower 
        third molar crown height; WM3: unworn lower third molar crown width) from Mendoza et al. (2002), measured in cm. Categorical variable (FS, i.e. 
        feeding style: B=browsers, G=grazers, M=mixed feeders) based on Pérez–Barbería and Gordon (2001). Phylogeny pruned from Hedges et al. (2015). 
        Taxonomic mismatches among these sources were resolved based on Wilson and Reeder (2005).        
        Hedges, S. B., J. Marin, M. Suleski, M. Paymer, and S. Kumar. 2015. Tree of life reveals clock-like speciation and diversification. 
        Molecular Biology and Evolution 32:835-845.        
        Mendoza, M., C. M. Janis, and P. Palmqvist. 2002. Characterizing complex craniodental patterns related to feeding behaviour in ungulates: 
        a multivariate approach. Journal of Zoology 258:223-246       
      Pérez–Barbería, F. J., and I. J. Gordon. 2001. Relationships between oral morphology and feeding style in the Ungulata: a phylogenetically
        controlled evaluation. Proceedings of the Royal Society of London. Series B: Biological Sciences 268:1023-1032.       
      Wilson, D. E., and D. M. Reeder. 2005. Mammal species of the world: A taxonomic and geographic reference. 
        Johns Hopkins University Press, Baltimore, Maryland.         
    
    4.2) Tree.tre 
        Ungulates' phylogeny, extracted from the mammalian phylogeny of         
        Hedges, S. B., J. Marin, M. Suleski, M. Paymer, and S. Kumar. 2015. Tree of life reveals clock–like speciation and diversification. Mol. Biol. Evol. 32:835–845.
        
    4.3) OUB.R, OUF.R, OUG.R
        R scripts for the analyses performed in the manuscript. Different files correspond to different regime setups of the feeding style variable.
        
    4.4) OU1.txt, OUB.txt, OUF.txt, OUG.txt
        Outputs of the model comparison conducted under the R scripts presented above (4.3). Different files correspond to different regime setups of the feeding style variable.
        
5) Ferula analyses

    In the models_ME directory there are input and output files from the mvSLOUCH analyzes of Ferula data with measurement error included, while in the models_no_ME                directory the analyzes of data without measurement error. In each directory, one can find the following files:

- input files: Data_ME.csv (with mesurment error) or Data_no_ME.csv (without measurement error) and tree file in Newick format (Ferula_fruits_tree.txt); the trait names in data files are abbreviated as follows: ration_canals – the proportion of oil ducts covering the space between median and lateral ribs, mean_gr_peri – periderm thickness, wings_area – wing area, thick_wings – wing thickness and seed_mass – seed mass,
- the results for 8 analyzed models (see Fig. 2 in the main text), each in separate directory named model1, model2 and so on,
- each model directory comprises the following files: two R scripts (for analyzes with diagonal and with upper triangular matrix Σyy; each model was run 1000 times), two csv files included information such as number of repetition (i), seed for preliminary analyzes generating starting point (seed_start_point), seed for the main analyses (seed) and AIC, AICc, SIC, BIC, R2 and loglik for each model run (these csv files are sorted according to AICc values), two directories containing results for 1000 analyzes, and two files extracted from these directories showing parameter estimation for the best models (with UpperTri and Diagonal matrix Σyy) 

6) test_rot_mvSLOUCH.R

    This R script demonstrates that the likelihood calculations of mvSLOUCH are rotation invariant.


    
