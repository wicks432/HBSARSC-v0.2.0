rm(list=ls())
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library(HBSARSC)

load("Simulation_1.rdata")

nbasis       = 20
iflagCenter  = 1
iflagHBsigma = 1
xmin         = 0
xmax         = 10
mcmc         = list(nblow0=1000,nblow=40000,smcmc=10000,nskip=1)

set.seed(1)
fout = hbsar(y=ydata, v=vdata, w=wdata, x=xdata, z=zdata, nbasis=nbasis, id_group=id_group, nint=500, 
             mcmc=mcmc, shape="Increasing", xmin=xmin, xmax=xmax,
             iflagCenter=iflagCenter, iflagHBsigma=iflagHBsigma)

ngroups   = fout$ngroup
xgrid     = fout$xgrid
f0xgridm  = fout$post.est$f0xgridm
f0Bxgridm = fout$post.est$f0Bxgridm
fxgridm   = fout$post.est$fxgridm
fxgridmm  = fout$post.est$fxgridmm
fxobsm    = fout$post.est$fxobsm
yresid    = ydata - fout$post.est$valpham - fout$post.est$wbetam

## Figure 2 (a)
matplot(xgrid,cbind(fxgridm,f0xgridm),type="l",xlab="x",ylab="f",
        lwd=c(rep(1,ngroups),3),lty=c(rep(5,ngroups),1),
        col=c(rep("black",ngroups),"red"))
legend("bottomright", legend=c("Upper Level","Lower Level"),col=c("red","black"),lty=c(1,5),lwd=c(4,2))

## Figure 2 (b)
matplot(xgrid,cbind(f0Bxgridm,f0xgridm,fxgridmm),xlab="X",ylab="Upper f",
        type="l",lwd=3,col=c("red","black","green"),lty=c(2,1,5))
legend("bottomright",legend=c("Posterior Mean","Biased Corrected","Average"),
       col=c("red","black","green"),lty=c(2,1,5),lwd=3)

## Figure 2 (c) & (d)
gid = c(18, 78)
ugroups = unique(fout$group)
for (j in gid) {

idx     = fout$group==j
yresidj = yresid[idx]
ydataj  = ydata[idx]
xdataj  = xdata[idx]
fxj     = fxgridm[,j]
f0j     = f0xgridm
fxobsj  = fxobsm[idx]
ymax    = ceiling(max(c(yresidj,fxj,fxgridt[,j],f0j,fxj)))
ymin    = floor(min(c(yresidj,fxj,fxgridt[,j],f0j,fxj)))
bout    = cbind(fxj,fxgridt[,j],f0j)

matplot(xgrid,bout,type="l",ylim=c(ymin,ymax),
        xlab="x",ylab="f",lwd=3,main=ugroups[j])
points(xdataj,yresidj,pch=19)
legend("bottomright",legend=c("Parametric Residuals","Lower-Level f","True f","Unbiased Upper-Level f"),
       col=c("black","black","red","green"),lty=c(0,1,2,3),lwd=c(0,2,2,2),pch=c(19,-1,-1,-1))

}
