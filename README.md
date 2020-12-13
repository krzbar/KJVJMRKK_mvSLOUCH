These are the R scripts and numerical results accompanying Bartoszek, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: Model comparison for multivariate Ornstein-Uhlenbeck-based models of trait evolution on large phylogenies".

The R setup for the manuscript was as follows: R version 3.6.1 (2019-09-12) Platform: x86_64-pc-linux-gnu (64-bit) Running under: openSUSE Leap 42.3

The exact output can depend on the random seed. However, in the script we have the option of rerunning the analyses as it was in the manuscript, i.e.
the random seeds that were used to generate the results are saved, included and can be read in.

The code is divided into several directories with scripts, random seeds and result files.

1) LikelihoodTesting

    Directory contains the script test_rotation_invariance_mvSLOUCH.R that demonstrates that mvSLOUCH's likelihood calculations are rotation invariant.
        
2) Carnivorans

    Directory contains files connected the the Carnivrons' vignette in mvSLOUCH
    
    2.1) Carnivora_mvSLOUCH_objects.RData    
           Full output of the running the R code in the vignette. With mvSLOUCH is a very bare-minimum subset of this file that allows for the creation of the            vignette.
    
3) SimulationStudy

    Directory contains all the output of the simulation study presented in the manuscript and scripts that allow for replication (the random number generator seeds are also provided) or running ones own simulation study, and scripts to generate graphs, and model comparison summary.
    
