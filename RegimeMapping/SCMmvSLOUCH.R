## This file accompanies the manuscript: 
## Bartoszek, Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## The code here generates the output for 
## Appendix SC7: Importing stochastic character mapping to mvSLOUCH: example analysis
## We show how to use tegime under stochastic character mapping for mvSLOUCH 

library(ape)
library(geiger)
library(PCMBaseCpp)
library(mvSLOUCH)
library(phytools)
library(rgl)

set.seed(95)


##################
### SIMULATION ###
##################

# Simulate tree
tree<-rtree(30)
png("S1.png")
plot(tree, font = 2, edge.width = 3, cex = 1.2)
add.scale.bar(length = 0.2, cex = 1, lwd = 2)
dev.off()
mvStree<-phyltree_paths(tree)

# Set up simulation
OUBMparameters<-list(vY0=matrix(c(5,2),ncol=1,nrow=2),
                     A=rbind(c(6,0),c(0,2)),
                     B=matrix(c(2,-1),ncol=1,nrow=2),
                     mPsi=matrix(c(5,2),ncol=1,nrow=2),
                     Syy=rbind(c(1,0.3),c(0,1)),vX0=matrix(0,1,1),
                     Sxx=matrix(1,1,1),Syx=matrix(0,ncol=1,nrow=2),
                     Sxy=matrix(0,ncol=2,nrow=1))
# Simulate data
OUBMdata<-simulMVSLOUCHProcPhylTree(mvStree,OUBMparameters)
# Check alignment
row.names(OUBMdata) == tree$tip.label

# Visualize continuous values on tips as colors
rbPal <- colorRampPalette(c('red','blue'))
Colors <- rbPal(30)[as.numeric(cut(OUBMdata[,"trait_3"],breaks = 30))]
png("S3a.png")

hist(OUBMdata[,"trait_3"],breaks = 20, main = NULL, xlab = "BM trait", 
     cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5,col="gray")
arrows(-1.25, 1, -1.25, 0, col = "red", lwd = 3)
dev.off()

## Discretization
## our aim is to create regimes based on the BM trait
Discr<-OUBMdata[,"trait_3"]
Discr<-ifelse(Discr < (-1),"B","A")


png("S3b.png")
plot(tree, label.offset = 0.2, font = 2, edge.width = 3, cex = 1.2)
tiplabels(pch = 21, bg = Colors, adj = 0.6, cex = 2)
# Add legend
Legend = rep(NA, 31)
Legend[c(1,16,31)] = c(as.numeric(gsub(".*,", "", gsub("[]]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[length(OUBMdata[,"trait_3"])]))),
                      median(c(as.numeric(gsub(",.*$", "", gsub("[(]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[1]))),as.numeric(gsub(".*,", "", gsub("[]]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[length(OUBMdata[,"trait_3"])]))))),
                      as.numeric(gsub(",.*$", "", gsub("[(]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[1]))))
legend("bottomleft", 
#x = -0.05, y = 6,
       legend = Legend,
       fill = colorRampPalette(colors = c('blue','red'))(31),
       border = NA,
       y.intersp = 0.1,
       cex = 1, text.font = 2)
dev.off()

# Graphical exploration of simulated data
tryCatch({
    if (rgl.useNULL()){
	print("There seems to be some problem with rgl, resorting to widgets.")
        options(rgl.printRglwidget = TRUE)
	RGLwidgetTraits3D<-rgl::plot3d(OUBMdata[,"trait_3"], OUBMdata[,"trait_1"], OUBMdata[,"trait_2"], 
        xlab = "Trait 3", ylab = "Trait 1", zlab = "Trait 2", 
	   size = 4,        cex = 4)
        htmltools::save_html(RGLwidgetTraits3D,file="S2.html")
    }else{
	rgl::plot3d(OUBMdata[,"trait_3"], OUBMdata[,"trait_1"], OUBMdata[,"trait_2"], 
    	    xlab = "Trait 3", ylab = "Trait 1", zlab = "Trait 2", size = 4,
    	    cex = 4)
        rgl.snapshot('S2.png', fmt = 'png')
    }
}, error=function(e){print(paste0("Cannot create 3D plot of traits: ",e))})

####################################
### STOCHASTIC CHARACTER MAPPING ###
####################################

name.check(tree,Discr)

# Model comparison #
# 1. Equal-rates
mtreeER <- make.simmap(tree, Discr, Q="empirical", model = "ER")
# 3. All-rates-different
mtreeARD <- make.simmap(tree, Discr, Q="empirical", model = "ARD")

# Likelihood ratio test
LR<-2*(mtreeARD$logL-mtreeER$logL)                                        # Likelihood ratio
p_value=1-pchisq(LR,df=attr(mtreeARD$logL,"df")-attr(mtreeER$logL,"df"))  # Significance


# Mapping
cols <- setNames(c("blue", "red"), sort(unique(Discr)))
png("VisualizeSimmap.png") ## plot not in SI
plotSimmap(mtreeER, cols, pts = FALSE, lwd = 3, fsize=0.8)
add.scale.bar(length=0.2,cex=0.7)
add.simmap.legend(colors = cols, vertical = TRUE, prompt = FALSE,
                  x = 4.5, y = 25)
dev.off()
# Description of SCM
describe.simmap(mtreeER)

# Simulate 500 stochastic character maps
mtrees <- make.simmap(tree, Discr, Q="empirical", model = "ER", nsim = 500)
XX <- describe.simmap(mtrees)

# Plot posterior probabilities on nodes of the tree
png("S5.png")
tree2<-tree
tree2$tip.label<-rep("",length(tree$tip.label))
plotTree(tree2, ftype = "b", fsize = 1.3, lwd = 3)
nodelabels(pie = XX$ace, piecol = cols, cex = 0.8)
add.scale.bar(length = 0.2, cex = 1, lwd = 2)
add.simmap.legend(colors = cols, vertical = TRUE, prompt = FALSE, 
                  x = 0.3, y = 5, fsize = 1.5)
tiplabels(pch = 21, bg = Colors, adj = 0.6, cex = 2)
# Add heatmap of BM tip trait
Legend = rep(NA, 31)
Legend[c(1,16,31)] = c(as.numeric(gsub(".*,", "", gsub("[]]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[length(OUBMdata[,"trait_3"])]))),
                      median(c(as.numeric(gsub(",.*$", "", gsub("[(]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[1]))),as.numeric(gsub(".*,", "", gsub("[]]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[length(OUBMdata[,"trait_3"])]))))),
                      as.numeric(gsub(",.*$", "", gsub("[(]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[1]))))
legend("bottomright",
       legend = Legend,
       fill = colorRampPalette(colors = c('blue','red'))(31),
       border = NA,
       y.intersp = 0.1,
       cex = 1, text.font = 2)

dev.off()
# Specify regime based on posterior probabilities of stochastic character mapping
#regimes<-rep("A",58);regimes[c(5:17,54:56)]<-"B"
regimes<-rep("A",nrow(tree.edge))
vBinternalnodes<-as.numeric(names(which(XX$ace[,"A"]<XX$ace[,"B"])))
vBinternalbranches<-which(tree$edge[,2]%in%vBinternalnodes)
vBtipnodes<-which(tree$tip.label%in%names(which(XX$tips[,"A"]<XX$tips[,"B"])))
vBtipbranches<-which(tree$edge[,2]%in%vBtipnodes)
vBbranches<-sort(c(vBinternalbranches,vBtipbranches))
regimes[vBbranches]<-"B"


# Confirm regime assignation
reg.col<-regimes
reg.col[reg.col=="A"]<-"blue"
reg.col[reg.col=="B"]<-"red"
png("S6.png")
plot(tree, cex = 1.2, font = 2, label.offset=0.2, edge.width=3, 
     edge.color = reg.col, show.tip.label=FALSE) 
tiplabels(pch = 21, bg = Colors, adj = 0.6, cex = 2)     
## below legends not needed when presented, in SI, side by side
#add.scale.bar(length=0.2, cex=1, lwd = 2);add.simmap.legend(colors = cols, vertical = T, prompt = FALSE, x = 0.3, y = 5, fsize = 1.5)
# Add heatmap of BM tip trait
##Legend = rep(NA, 31);Legend[c(1,16,31)] = c(as.numeric(gsub(".*,", "", gsub("[]]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[length(OUBMdata[,"trait_3"])]))),median(c(as.numeric(gsub(",.*$", "", gsub("[(]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[1]))),as.numeric(gsub(".*,", "", gsub("[]]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[length(OUBMdata[,"trait_3"])]))))),as.numeric(gsub(",.*$", "", gsub("[(]", "", levels(cut(OUBMdata[,"trait_3"],breaks = 30))[1]))));legend("bottomright",legend = Legend,fill = colorRampPalette(colors = c('blue','red'))(31),border = NA,y.intersp = 0.1,cex = 1, text.font = 2)
dev.off()


################
### ANALYSIS ###
################


# Preparing the subdataset
mvData<-data.matrix(OUBMdata[,c("trait_1","trait_2")])

# Model comparison #
# Brownian motion
BMestim<-BrownianMotionModel(mvStree,mvData)

# Global OU 
OU1estim<-ouchModel(mvStree,mvData,regimes = NULL,
                    Atype = "Diagonal", Syytype = "UpperTri",
                    diagA = "Positive")

# Regime OU
OU2estim<-ouchModel(mvStree, mvData, regimes = regimes, root.regime = "A",
                    Atype = "Diagonal", Syytype = "UpperTri",
                    diagA = "Positive")

###############
### OUTPUTS ###
###############

cat('Transition rates (stochastic character mapping):','\n',
    '  Equal:','\n',mtreeER$Q,'\n',
    '  Unequal:','\n',mtreeARD$Q,'\n','\n',
    'Likelihood ratio test:','\n',
    '  LR =',LR[1],'\n',
    '  df =',attr(mtreeARD$logL,"df")-attr(mtreeER$logL,"df"),'\n',
    '  p value =',p_value[1],'\n','\n',
    'Model comparison under mvSLOUCH (AICc):','\n',
    '  BM =',BMestim$ParamSummary$aic.c,'\n',
    '  OU1 =',OU1estim$FinalFound$ParamSummary$aic.c,'\n',
    '  OU2 =',OU2estim$FinalFound$ParamSummary$aic.c,'\n')
