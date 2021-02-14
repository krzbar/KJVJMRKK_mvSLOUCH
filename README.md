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
        Phenotypic measurements of muzzle width (MZW, measured at the outer junction of the boundary between the maxilla and premaxilla), unworn lower third molar crown height (HM3 , measured from the base of the crown to the tip of the protoconid), width (WM3 , measured at the occlusal surface, between the outer aspects of the protoconid and the entoconid), and feeding style (FS, three levels B = Browser; M = Mixed feeder; G = Grazer) for 104 Ungulate species. Measurements are compiled from  
        Mendoza, M., C. M. Janis, and P. Palmqvist. 2002. Characterizing complex craniodental patterns related to feeding behaviour in ungulates: a multivariate approach. J. Zool. 258:223–246.        
        Pérez–Barberı́a, F. J. and I. J. Gordon. 2001. Relationships between oral morphology and feeding style in the ungulata: a phylogenetically controlled evaluation. P. Roy. Soc. Lond. B Bio. 268:1023–1032.
    
    4.2) OUB.R, OUF.R, OUG.R
        R scripts for the analyses performed in the manuscript. Different files correspond to different regime setups of the feeding style variable.
        
    4.3) Tree.tre 
        Ungulates' phylogeny, extracted from the mammalian phylogeny of         
        Hedges, S. B., J. Marin, M. Suleski, M. Paymer, and S. Kumar. 2015. Tree of life reveals clock–like speciation and diversification. Mol. Biol. Evol. 32:835–845.

    
