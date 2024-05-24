## This file accompanies the manuscript: 
## Bartoszek, Clarke, Fuentes Gonzalez, Mitov, Pienaar, Piwczynski, Puchalka, Spalik and Voje " Fast mvSLOUCH: multivariate Ornstein–Uhlenbeck-based models of trait evolution on large phylogenies"
## It generates Fig. 1 of the Supplementary Material.

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## modified solution of Owlright from
##https://stackoverflow.com/questions/14196551/add-isoclines-and-or-direction-field-to-plot

library(deSolve) 
library(phaseR)
library(latex2exp)


## ======================================================================
## edit these two variables to control the aesthetics of the flow field
ff_points<-15 ## controls the density of points for the flow field
ff_frac<-1.1 ## controls length of arrows for the flow field
## ======================================================================

f_titlematrix<-function(lA){
 #TeX(r'($\left[\begin{array}{cc} 1 & 1 \\ 1 & 1 \end{array} \right]$)')
 paste0(eval(expression(signif(lA$A[,1],3)))," ",eval(expression(signif(lA$A[,2],3))))
}




model.OU <- function(t, y, parameters){
  with(as.list(parameters),{

    OU1<-y[1] 
    OU2<-y[2]
    dOU1 <- (-1)*( a11*(OU1-theta1) + a12*(OU2-theta2))
    dOU2 <- (-1)*( a21*(OU1-theta1) + a22*(OU2-theta2))

    list(c(dOU1,dOU2))
  })
}

lA<-list(
list(A=NA,eigA=list(values=c(-0.01,-0.05),vectors=cbind(c(1,1),c(-1,5)))),
list(A=NA,eigA=list(values=c(0.01,0.05),vectors=cbind(c(1,1),c(-1,5)))),
list(A=NA,eigA=list(values=c(0.01,-0.05),vectors=cbind(c(1,1),c(-1,5)))),
list(A=NA,eigA=list(values=c(-0.01,0.05),vectors=cbind(c(1,1),c(-1,5)))),
#list(A=NA,eigA=list(values=c(0.01,0.01),vectors=cbind(c(1,1),c(-1,5)))),
#list(A=NA,eigA=list(values=c(-0.01,-0.01),vectors=cbind(c(1,1),c(-1,5)))),
#list(A=NA,eigA=list(values=c(0.01,-0.01),vectors=cbind(c(1,1),c(-1,5)))),
#list(A=NA,eigA=list(values=c(-0.01,0.01),vectors=cbind(c(1,1),c(-1,5)))),
#list(A=NA,eigA=list(values=c(-0.01,-0.05),vectors=cbind(c(0,1),c(-1,0)))),
#list(A=NA,eigA=list(values=c(0.01,0.05),vectors=cbind(c(0,1),c(-1,0)))),
#list(A=NA,eigA=list(values=c(0.01,-0.05),vectors=cbind(c(0,1),c(-1,0)))),
#list(A=NA,eigA=list(values=c(-0.01,0.05),vectors=cbind(c(0,1),c(-1,0)))),
list(A=NA,eigA=list(values=c(0.01,0.01),vectors=cbind(c(0,1),c(-1,0)))),
list(A=NA,eigA=list(values=c(-0.01,-0.01),vectors=cbind(c(0,1),c(-1,0)))),
list(A=NA,eigA=list(values=c(0.01,-0.01),vectors=cbind(c(0,1),c(-1,0)))),
list(A=NA,eigA=list(values=c(-0.01,0.01),vectors=cbind(c(0,1),c(-1,0)))),
list(A=cbind(c(0.01,0),c(0.02,0.05)),eigA=NA),
list(A=cbind(c(0.05,0),c(0.02,0.01)),eigA=NA),
list(A=cbind(c(0.01,0),c(0.005,0.05)),eigA=NA),
list(A=cbind(c(0.05,0),c(0.005,0.01)),eigA=NA),
list(A=cbind(c(0.01,0.02),c(0,0.05)),eigA=NA),
list(A=cbind(c(0.05,0.02),c(0,0.01)),eigA=NA),
list(A=cbind(c(0.01,0.005),c(0,0.05)),eigA=NA),
list(A=cbind(c(0.05,0.005),c(0,0.01)),eigA=NA),
list(A=cbind(c(0.01,0),c(0.02,0.01)),eigA=NA),
list(A=cbind(c(0.01,0),c(0.005,0.01)),eigA=NA),
list(A=cbind(c(0.01,0.02),c(0,0.01)),eigA=NA),
list(A=cbind(c(0.01,0.005),c(0,0.01)),eigA=NA)
)

lA<-sapply(lA,function(x){
    if (is.na(x$A[1])){
	P<-x$eigA$vectors
	L<-diag(x$eigA$values)
	x$A<-P%*%L%*%solve(P)
    }else{
	x$eigA<-eigen(x$A)	
    }
    x
},simplify=FALSE)

vi_toplot<-c(
2,1,5,6,
9,17,10,12,
15,20,16,14)
fig_size<-240*1.125


m_y0<-rbind(
c(-0.025,0.025,0.025,-0.025), #1 (1,2)
c(-0.05,0.05,0.05,-0.05), #2 (1,1)
c(-1,1,1,-1), #3 X
c(-1,1,1,-1), #4 X
c(-0.05,0.05,0.05,-0.05), #5 (1,3)
c(-0.05,0.05,0.05,-0.05), #6 (1,4)
c(-1,1,1,-1), #7 X
c(-0.1,0.1,0.1,-0.1), #8 X
c(-0.075,0.075,0.075,-0.075), #9 (2,1)
c(-0.075,-0.075,0.075,0.075), #10 (2,3)
c(-1,1,1,-1), #11 X
c(-0.075,0.075,0.075,-0.075), #12 (2,4)
c(-1,1,1,-1), #13 X
c(-0.075,0.075,0.075,-0.075), #14 (3,4)
c(-0.075,0.075,0.075,-0.075), #15 (3,1)
c(-0.075,0.075,0.075,-0.075), #16 (3,3)
c(-0.075,0.075,0.075,-0.075), #17 (2,2)
c(-0.1,0.1,0.1,-0.1), #18 X
c(-1,1,1,-1), #19 X
c(-0.075,-0.075,0.075,0.075) #20 (3,2)
)

png("Fig1_expdynamics_traj.png",width = 4*fig_size, height = 3*fig_size )
par(mfrow=c(3,4))
sapply(1:length(vi_toplot),function(j,lA,vi_toplot,m_y0){
    i<-vi_toplot[j]
    mA<-lA[[i]]$A
    params.OU<-c(a11=mA[1,1],a12=mA[1,2],a21=mA[2,1],a22=mA[2,2],theta1=0,theta2=0)
    data.OU<-as.data.frame(lsoda(c(OU1=0,OU2=0),seq(1,250,by=0.01), model.OU, params.OU))
    # plot the trajectories of the system
    plot(data.OU$OU1, data.OU$OU2, type="l", col="blue", main=NULL, # main=f_titlematrix(lA[[i]]),
        xlab="", ylab="", xlim=c(-0.1,0.1), ylim=c(-0.1,0.1),cex.axis=1.5,font.axis=2,font.lab=2,font.main=2,cex.main=2)
    title(paste0(letters[j],")"),adj=0,font.main=2,cex.main=2)
    #add Nullclines
    nullclines(model.OU, xlim=c(-0.1,0.1),ylim=c(-0.1,0.1), parameters=params.OU, system="two.dim", col=c("blue","blue"), add=TRUE,add.legend=FALSE,lwd=2)
    flowField(model.OU, xlim=c(-0.1,0.1),ylim=c(-0.1,0.1), parameters=params.OU, system="two.dim", add=TRUE,add.legend=FALSE,frac=ff_frac,points=ff_points) #, col=c("black")
    trajectory(model.OU,tlim=c(0,1000),y0=rbind(m_y0[i,1:2],m_y0[i,3:4]), parameters=params.OU, system="two.dim",frac=ff_frac)
},lA=lA,vi_toplot=vi_toplot,m_y0=m_y0,simplify=FALSE)
dev.off()


