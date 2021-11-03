C=12.5

res_lin=readRDS(file=paste("mcmc/monthly/MCMC_linear_20k_",C,".rds",sep=""))
res_3=readRDS(file=paste("mcmc/monthly/MCMC_RCS3_20k_",C,".rds",sep=""))
res_7=readRDS(file=paste("mcmc/monthly/MCMC_RCS7_20k_",C,".rds",sep=""))
res_LGLL=readRDS(file="mcmc/monthly/MCMC_LGLL_20k.rds")

library(readxl)
library(splines)
library(RcppStorage)
library(mvtnorm)
dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])
datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])
##################################################################
##################################################################
##################################################################
npar=10000
nburnin=2000
nMCMC = 20000+nburnin

datnames=c("natgas","coffee","cotton","aluminium")
n=length(datnames)

library(LaplacesDemon)
options(scipen=999)
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
a=0
r = ((1.05)^(1/12))-1
v_x=1

c_s0 = 0.01
c_df = 10
c_alpha = 2
c_beta = 20
c_detvar=20^2

ngrid = 20
nint = 128
int.grid = seq(from=-4.0,to=4.0,length.out=nint)
int.wts = dnorm(int.grid)*0.5*(int.grid[2]-int.grid[1])
xleft=0.0
xright=C
grid=seq(from=0,to=C,length.out=ngrid)
vals=seq(from=0,to=C,length.out=ngrid)
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
logpost=function(chainres,postmean,ll,FUN)
{
  theta=chainres$theta[(nburnin+1):nMCMC,]
  loglike=chainres$loglike[(nburnin+1):nMCMC]
  qCov = chainres$q_cov
  
  M=nMCMC-nburnin
  numerator_log=numeric(M)
  denominator_log=numeric(M)
  newlog=numeric(M)
  
  for(i in 1:M)
  {
    numerator_log[i] = min(0,ll-loglike[i])+dmvnorm(postmean,mean=theta[i,],sigma=qCov,log=TRUE)
    theta_j = rmvnorm(1,mean=postmean,sigma=qCov)
    newlog[i]=FUN(theta_j)
    denominator_log[i]=min(0,newlog[i]-ll)
  }
  
  num_max = max(numerator_log)
  den_max = max(denominator_log)
  
  tmp1=num_max+log(mean(exp(numerator_log-num_max)))
  tmp2=den_max+log(mean(exp(denominator_log-den_max)))
  
  return(list(nl=numerator_log,dl=denominator_log,lp=tmp1-tmp2,newlog=newlog))
}
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
logl_stoch = function(pars)
{
  ll=SMC_sample(pars,p,npar=10000,N,xleft,xright,grid,vals,a,C,r,v_x=1,int.grid,int.wts,nint,ngrid,30,0.5)$loglike_noprior
  
  v_prior = -(c_df*pars[1])+0.5*c_df*log(0.5*c_s0*c_df)-log(gamma(0.5*c_df))+log(2)-0.5*(c_df*c_s0*exp(-2.0*pars[1]))
  delta=0.5+0.5*tanh(pars[2])
  delta_prior = -log(gamma(c_alpha))-log(gamma(c_beta))+log(gamma(c_alpha+c_beta))+log(delta)*(c_alpha-1.0)+log(1.0-delta)*(c_beta-1.0)+log(abs(-0.5+0.5*(2.0*delta-1.0)^2))
  b_prior = -0.5*log(2*pi)-0.5*(pars[3])^2
  
  return(ll+v_prior+delta_prior+b_prior)
}
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
logl_linear <- function(kpars)
{
  np=length(kpars)
  delta=0.5+0.5*tanh(kpars[np-1])
  b=exp(kpars[np])
  beta = (1.0-delta)/(1.0+r)
  
  rep= Storagefunc(xleft,xright,grid,vals,a,b, C, delta, beta,int.grid,int.wts,nint,ngrid, 30)
  
  x = numeric(N)
  x[1] = finv(exp(p[1]-k(1,kpars[1:(np-2)])),rep$xleft,rep$xright,rep$grid,rep$vals,a,b,C)
  
  loglike=0
  
  for(j in 2:N)
  {
    tmp=exp(p[j]-k(j,kpars[1:(np-2)]))
    x[j]=finv(tmp,rep$xleft,rep$xright,rep$grid,rep$vals,a,b,C)
    z =x[j]-(1-delta)*evalSigma(x[j-1],rep$xleft,rep$xright,rep$grid,rep$vals,C)
    Jacobian = abs(tmp/(evalf_dx(x[j], rep$xleft, rep$xright, rep$grid, rep$vals, a,b,C)))
    loglike = loglike + log((1/sqrt(2*pi))*exp(-0.5*z^2)*Jacobian)
  }
  
  delta_prior = -log(gamma(c_alpha))-log(gamma(c_beta))+log(gamma(c_alpha+c_beta))+log(delta)*(c_alpha-1.0)+log(1.0-delta)*(c_beta-1.0)+log(abs(-0.5+0.5*(2.0*delta-1.0)^2))
  b_prior = -0.5*log(2*pi)-0.5*(kpars[np])^2
  det_prior= sum(-0.5*log(2*pi*c_detvar)-((kpars[1:(np-2)])^2)/(2*c_detvar))
  
  return(loglike+delta_prior+b_prior+det_prior)
}

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
logl_RCS <- function(kpars)
{
  np=length(kpars)
  delta=0.5+0.5*tanh(kpars[np-1])
  b=exp(kpars[np])
  beta = (1.0-delta)/(1.0+r)
  
  rep= Storagefunc(xleft,xright,grid,vals,a,b, C, delta, beta,int.grid,int.wts,nint,ngrid, 30)
  
  x = numeric(N)
  x[1] = finv(exp(p[1]-k(1,kpars[1:(np-2)])),rep$xleft,rep$xright,rep$grid,rep$vals,a,b,C)
  
  loglike=0
  
  for(j in 2:N)
  {
    tmp=exp(p[j]-k(j,kpars[1:(np-2)]))
    x[j]=finv(tmp,rep$xleft,rep$xright,rep$grid,rep$vals,a,b,C)
    z =x[j]-(1-delta)*evalSigma(x[j-1],rep$xleft,rep$xright,rep$grid,rep$vals,C)
    Jacobian = abs(tmp/(evalf_dx(x[j], rep$xleft, rep$xright, rep$grid, rep$vals, a,b,C)))
    loglike = loglike + log((1/sqrt(2*pi))*exp(-0.5*z^2)*Jacobian)
  }
  
  delta_prior = -log(gamma(c_alpha))-log(gamma(c_beta))+log(gamma(c_alpha+c_beta))+log(delta)*(c_alpha-1.0)+log(1.0-delta)*(c_beta-1.0)+log(abs(-0.5+0.5*(2.0*delta-1.0)^2))
  b_prior = -0.5*log(2*pi)-0.5*(kpars[np])^2
  det_prior= sum(-0.5*log(2*pi*c_detvar)-((kpars[1:(np-2)])^2)/(2*c_detvar))
  
  return(loglike+delta_prior+b_prior+det_prior)
}
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
kalman.filter = function(sysq,sxsq,phi,mu1,omegasq,y)
{
  T = length(y)
  predmean = vector(length=T,mode="numeric")
  predvar = predmean
  filtmean = predmean
  filtvar = predmean
  loglikelihood = predmean
  resid = predmean
  
  # first time step
  predmean[1] = mu1
  predvar[1] = omegasq
  
  filtmean[1] = (predmean[1]*sysq + y[1]*predvar[1])/(sysq + predvar[1])
  filtvar[1] = sysq*predvar[1]/(sysq+predvar[1])
  
  loglikelihood[1] = dnorm(x=y[1],mean=predmean[1],sd=sqrt(predvar[1]+sysq),log=TRUE)
  resid[1] = (y[1]-predmean[1])/sqrt(predvar[1]+sysq)
  
  # the remaining time steps
  for(t in 2:T)
  {
    predmean[t] = phi*filtmean[t-1]
    predvar[t] = phi^2*filtvar[t-1]+sxsq
    
    filtmean[t] = (predmean[t]*sysq + y[t]*predvar[t])/(sysq + predvar[t])
    filtvar[t] = sysq*predvar[t]/(sysq+predvar[t])
    
    loglikelihood[t] = dnorm(x=y[t],mean=predmean[t],sd=sqrt(predvar[t]+sysq),log=TRUE)
    resid[t] = (y[t]-predmean[t])/sqrt(predvar[t]+sysq)
  }
  # return generated quantities
  return(list(predmean=predmean,predvar=predvar,filtmean=filtmean,
              filtvar=filtvar,loglikelihood=loglikelihood,resid=resid))
}

logl_LGLL = function(pars)
{
  sysq = exp(2*pars[1]) #log(b) -> b^2      
  sxsq = exp(2*pars[2]) #log(v) -> v^2
  ret = kalman.filter(sysq=sysq,sxsq=sxsq,phi=1.0,mu1=0.0,omegasq = 1000000.0,y=p)
  v_prior = -(c_df*pars[2])+0.5*c_df*log(0.5*c_s0*c_df)-log(gamma(0.5*c_df))+log(2)-0.5*(c_df*c_s0*exp(-2.0*pars[2]))
  b_prior = -0.5*log(2*pi)-0.5*(pars[1])^2
  return(sum(ret$loglikelihood[2:length(p)])+v_prior+b_prior)
}

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# STOCHASTIC

for(j in 1:n)
{
  p = log(datlist[[j]])
  N=length(p)
  #############################
  #############################
  #############################
  # Stochastic
  res_stoch=readRDS(file=paste("mcmc/monthly/",datnames[j],"_",10000,"_",C,".rds",sep=""))
  stoch_postmean=colMeans(res_stoch$theta[(nburnin+1):nMCMC,])
  
  ll = logl_stoch(stoch_postmean)
  
  lp=logpost(res_stoch,stoch_postmean,ll,logl_stoch)
  marglike=list(marglike=ll-lp$lp,lp=lp,ll=ll)
  
  
  saveRDS(marglike, file=paste("mcmc/monthly/margloglike_stoch_",datnames[j],"_",C,".rds",sep=""))
}


##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# LINEAR
marglike = list()

for(j in 1:n)
{
  print(j)
  p = log(datlist[[j]])
  N=length(p)
  #############################
  #############################
  #############################
  # Linear
  k = function(t,kpars)
  {
    return(kpars[1]+kpars[2]*(t/N))
  }
  
  t_lin=res_lin[[j]]$theta[(nburnin+1):nMCMC,]
  lin_postmean=colMeans(t_lin)

  ll=logl_linear(lin_postmean)
  
  lp=logpost(res_lin[[j]],lin_postmean,ll,logl_linear)
  marglike[[j]]=list(marglike=ll-lp$lp,lp=lp,ll=ll)
}

saveRDS(marglike, file=paste("mcmc/monthly/margloglike_linear_",C,".rds",sep=""))

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# RCS3
marglike = list()

for(j in 1:n)
{
  print(j)
  p = log(datlist[[j]])
  N=length(p)
  #############################
  #############################
  #############################
  # RCS3
  knots=ns((1:N)/N, knots=quantile((1:N)/N,probs=c(0.25,0.5,0.75)),intercept = TRUE)
  
  k = function(t,kpars)
  {
    return(sum(knots[t,]*kpars))
  }
  
  t_3=res_3[[j]]$theta[(nburnin+1):nMCMC,]
  postmean_3=colMeans(t_3)
  
  ll=logl_RCS(postmean_3)
  
  lp=logpost(res_3[[j]],postmean_3,ll,logl_RCS)
  marglike[[j]]=list(marglike=ll-lp$lp,lp=lp,ll=ll)
}

saveRDS(marglike, file=paste("mcmc/monthly/margloglike_RCS3_",C,".rds",sep=""))

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# RCS7
marglike = list()

for(j in 1:n)
{
  print(j)
  p = log(datlist[[j]])
  N=length(p)  
  #############################
  #############################
  #############################
  # RCS7
  knots=ns((1:N)/N, knots=quantile((1:N)/N,probs=c(0.125,0.250,0.375,0.500,0.625,0.750,0.875)),intercept = TRUE)
  
  k = function(t,kpars)
  {
    return(sum(knots[t,]*kpars))
  }
  
  t_7=res_7[[j]]$theta[(nburnin+1):nMCMC,]
  postmean_7=colMeans(t_7)
  
  ll=logl_RCS(postmean_7)
  
  lp=logpost(res_7[[j]],postmean_7,ll,logl_RCS)
  marglike[[j]]=list(marglike=ll-lp$lp,lp=lp,ll=ll)
  
}

saveRDS(marglike, file=paste("mcmc/monthly/margloglike_RCS7_",C,".rds",sep=""))

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# LGLL
marglike = list()

for(j in 1:n)
{
  p = log(datlist[[j]])
  N=length(p)
  #############################
  #############################
  #############################
  # LGLL
  t_LGLL=res_LGLL[[j]]$theta[(nburnin+1):nMCMC,]
  LGLL_postmean=colMeans(t_LGLL)
  ll=logl_LGLL(LGLL_postmean)
  
  lp=logpost(res_LGLL[[j]],LGLL_postmean,ll,logl_LGLL)
  marglike[[j]]=list(marglike=ll-lp$lp,lp=lp,ll=ll)
  
}

saveRDS(marglike, file=paste("mcmc/monthly/margloglike_LGLL_",C,".rds",sep=""))