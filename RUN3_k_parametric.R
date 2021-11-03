C=12.5

library(readxl)
library(splines)
dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])
datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])

lin=readRDS(file=paste("mcmc/monthly/MCMC_linear_20k_",C,".rds",sep=""))
RCS3=readRDS(file=paste("mcmc/monthly/MCMC_RCS3_20k_",C,".rds",sep=""))
RCS7=readRDS(file=paste("mcmc/monthly/MCMC_RCS7_20k_",C,".rds",sep=""))
##################################################################################
##################################################################################
##################################################################################
nburnin=2000
nMCMC = 20000+nburnin

klist=list()

for(i in 1:length(datlist))
{
  p = log(datlist[[i]])
  N=length(p)
  kmat=matrix(0,nrow=N,ncol=3)
  ########################################
  ########################################
  ########################################
  # Linear
  k = function(t,kpars)
  {
    return(kpars[1]+kpars[2]*(t/N))
  }
  linpars=colMeans(lin[[i]]$theta[(nburnin+1):nMCMC,])
  for(q in 1:N)
  {
    kmat[q,1]=k(q,linpars[1:2])
  }
  ########################################
  ########################################
  ########################################
  # RCS3
  knots=ns((1:N)/N, knots=quantile((1:N)/N,probs=c(0.25,0.5,0.75)),intercept = TRUE)
  k = function(t,kpars)
  {
    return(sum(knots[t,]*kpars))
  }
  RCS3pars=colMeans(RCS3[[i]]$theta[(nburnin+1):nMCMC,])
  for(q in 1:N)
  {
    kmat[q,2]=k(q,RCS3pars[1:5])
  }
  ########################################
  ########################################
  ########################################
  # RCS7
  knots=ns((1:N)/N, knots=quantile((1:N)/N,probs=c(0.125,0.250,0.375,0.500,0.625,0.750,0.875)),intercept = TRUE)
  k = function(t,kpars)
  {
    return(sum(knots[t,]*kpars))
  }
  RCS7pars=colMeans(RCS7[[i]]$theta[(nburnin+1):nMCMC,])
  for(q in 1:N)
  {
    kmat[q,3]=k(q,RCS7pars[1:9])
  }
  ########################################
  ########################################
  ########################################
  ########################################
  
  klist[[i]]=kmat
}

saveRDS(klist,file=paste("mcmc/monthly/kvals_parametric_",C,".rds",sep=""))

##################################################################################
##################################################################################
##################################################################################
