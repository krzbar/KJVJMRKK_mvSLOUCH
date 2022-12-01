## This file accompanies the manuscript: 
## Bartoszek, Tredgett Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczyński, Puchałka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

# ~~~~~~~~~~~~~~~~~~ Plot mPsi with regression CIs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DFs1 <- data.frame(
  Ecology=factor(colnames(best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)),
  Trait= rep( "plant height(m)", length(colnames(best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)) ),
  Estimated.Point=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point["plant.height",],
  upper=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Upper.end["plant.height",],
  lower=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Lower.end["plant.height",]
)

DFs2 <- data.frame(
  Ecology=factor(colnames(best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)),
  Trait= rep( "seed mass(mg)", length(colnames(best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)) ),
  Estimated.Point=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point["seed.mass",],
  upper=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Upper.end["seed.mass",],
  lower=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Lower.end["seed.mass",]
)

DFs3 <- data.frame(
  Ecology=factor(colnames(best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)),
  Trait= rep( "leaf area(mm2)", length(colnames(best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)) ),
  Estimated.Point=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point["leaf.area",],
  upper=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Upper.end["leaf.area",],
  lower=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Lower.end["leaf.area",]
)

DFs4 <- data.frame(
  Ecology=factor(colnames(best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)),
  Trait= rep( "leaf mass(mg)", length(colnames(best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point)) ),
  Estimated.Point=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Estimated.Point["leaf.mass",],
  upper=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Upper.end["leaf.mass",],
  lower=best_model$best_output$ParamSummary$confidence.interval$regression.summary$mPsi.regression.confidence.interval$Lower.end["leaf.mass",]
)


DFs <- rbind(DFs1, DFs2, DFs3, DFs4)
DFs$Ecology2 <-  as.factor(sapply(DFs$Ecology,function(x){as.numeric(strsplit(as.character(x),"_")[[1]][2])},simplify=TRUE))


ggplot_optima_by_Ellenberg <- ggplot( DFs, aes(x=Ecology2, y=Estimated.Point, group=Trait, color=Trait) ) + 
  geom_point(size=2.5) + 
  geom_errorbar(ggplot2::aes(ymin=lower, ymax=upper), width=0.2) + 
  theme_bw() + 
  theme( axis.line = element_line(size=0.1) ) + 
  ## X axis related 
  theme( axis.title.x = element_blank() ) + 
  theme( axis.ticks.x = element_line(size=0.4) ) + 
  theme( axis.text.x = element_text(size=12,face="bold") ) + 
  theme( axis.title.x = element_text(size=12, face="bold",margin = margin(t=8, r=0, b=0, l=0)) ) + 
  xlab( paste0("Ellenberg values") ) + 
  scale_x_discrete(position = "bottom")  + 
  ## Y axis related 
  theme( axis.ticks.y = element_line(size=0.4) ) + 
  theme( axis.text.y = element_text(size=12,face="bold") ) + 
  theme( axis.title.y = element_text(size=12,face="bold", margin = margin(t=0, r=8, b=0, l=0)) ) + 
  ylab( paste0("log(optimum value)") ) +  
  scale_y_continuous(breaks = round( seq( -4, 9, by=1), 1) ) + 
  
  ## Strip 
  theme( panel.border = element_blank() ) + 
  theme( panel.grid.major = element_blank() ) + 
  theme( panel.grid.minor = element_blank() ) + 
  theme( panel.grid.minor.x = element_blank() ) + 
  theme( panel.grid.minor.y = element_blank() ) + 
  theme( panel.ontop = element_blank() ) + 
  theme( plot.background = element_blank() ) + 
  
  ## Legend 
  theme( legend.title = element_blank() ) + 
  theme( legend.text = element_text(size=12,face="bold") ) + 
  theme( legend.spacing.x = unit(0.2, 'cm') ) + 
  theme( legend.position = "bottom", legend.box = "horizontal", legend.direction = "horizontal" ) + 
  theme( legend.margin = margin(t=-5, r=0, b=0, l=0)) +
  guides(fill = guide_legend(nrow = 1)) 

  ggsave( paste0("./",c_resultsdirectory,"/Optima_bestmodel.pdf"), 
          ggplot_optima_by_Ellenberg, width=190, height=150, units="mm", useDingbats=FALSE ) 

# ~~~~~~~~~~~~~~~~~~ End plot mPsi with regression CIs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
