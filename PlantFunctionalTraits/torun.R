## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .



source("usercontrolledvariables.R")
source("tosource.R")
source("datasetup.R")
source("regimes.R")
source("startingmodels.R")
source("model1.R")
source("model2.R")
source("model3.R")
source("model4.R")
source("model5.R")
source("model6.R")
source("do_parametric_bootstrap.R")
if (b_make_optima_plot){source("make_optima_plot.R")}
source("modelBM.R") ## the results on this dataset under BM have to be compared manually, they are not in the automatic comparison as they are substantially worse
