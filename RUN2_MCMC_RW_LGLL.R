library(readxl)
dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])
datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])

c_s0 = 0.01
c_df = 10
##########################################################################################
##########################################################################################
##########################################################################################
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

logl_RW = function(pars)
{
  v = exp(pars[1])
  v_prior = -(c_df*pars[1])+0.5*c_df*log(0.5*c_s0*c_df)-log(gamma(0.5*c_df))+log(2)-0.5*(c_df*c_s0*exp(-2.0*pars[1]))
  return(sum(dnorm(diff(p),mean=0.0,sd=v,log=TRUE))+v_prior)
}
##########################################################################################
##########################################################################################
##########################################################################################
MCMC = function(kpars,p,logl)
{
  npar=length(kpars)
  theta=matrix(0,nrow=Nmcmc,ncol=npar)
  theta[1,] = kpars
  
  loglike=numeric(Nmcmc)
  acc = numeric(Nmcmc-1)
  accprobs = numeric(Nmcmc-1)
  proposals = list()
  
  loglike[1]=logl(theta[1,])
  
  C_0 = diag(npar)*rep(0.001,npar)
  
  library(mvtnorm)
  
  for(t in 2:(Nburnin/2))
  {
    p_theta = rmvnorm(1,mean=theta[t-1,],sigma=C_0)
    p_loglike = logl(p_theta)
    
    if(!is.finite(p_loglike))
    {
      p_loglike=-5000;
    }
    
    accp = exp(min(0,p_loglike-loglike[t-1]))
    
    accprobs[t-1]=accp
    
    if(runif(1)<accp)
    {
      theta[t,]=p_theta
      loglike[t]=p_loglike
      acc[t-1]=1
    }else
    {
      theta[t,]=theta[t-1,]
      loglike[t]=loglike[t-1]
      accprobs[t-1]=0
    }
  }
  
  sd=(2.4^2)/npar
  
  for(t in ((Nburnin/2)+1):Nburnin)
  {
    if(npar==1)
    {
      C_t = sd*var(theta[Nburnin_adapt:(t-1),])+sd*0.000001*diag(npar)
    }else
    {
      C_t = sd*cov(theta[Nburnin_adapt:(t-1),])+sd*0.000001*diag(npar)
    }
    p_theta = rmvnorm(1,mean=theta[t-1,],sigma=C_t)
    p_loglike = logl(p_theta)
    
    if(!is.finite(p_loglike))
    {
      p_loglike=-5000;
    }
    
    accp = exp(min(0,p_loglike-loglike[t-1]))
    
    accprobs[t-1]=accp
    
    if(runif(1)<accp)
    {
      theta[t,]=p_theta
      loglike[t]=p_loglike
      acc[t-1]=1
    }else
    {
      theta[t,]=theta[t-1,]
      loglike[t]=loglike[t-1]
      accprobs[t-1]=0
    }
  }
  
  for(t in (Nburnin+1):Nmcmc)
  {
    p_theta = rmvnorm(1,mean=theta[t-1,],sigma=C_t)
    
    proposals[[t-1]]=list(theta[t-1,],C_t)
    p_loglike = logl(p_theta)
    
    if(!is.finite(p_loglike))
    {
      p_loglike=-5000;
    }
    
    accp = exp(min(0,p_loglike-loglike[t-1]))
    
    accprobs[t-1]=accp
    
    if(runif(1)<accp)
    {
      theta[t,]=p_theta
      loglike[t]=p_loglike
      acc[t-1]=1
    }else
    {
      theta[t,]=theta[t-1,]
      loglike[t]=loglike[t-1]
      accprobs[t-1]=0
    }
  }
  
  res = list(loglike=loglike,theta=theta,acc=acc,accprobs=accprobs,q_cov=C_t)
  
  return(res)
}
##########################################################################################
##########################################################################################
##########################################################################################
Nburnin=2000
Nburnin_adapt=500
Nmcmc = 20000+Nburnin
##########################################################################################
##########################################################################################
##########################################################################################
res_RW=list()
res_LGLL=list()

for(i in 1:4)
{
  p = log(datlist[[i]])
  N=length(p)

  #RW
  RW_init=sd(diff(p)) #v
  res_RW[[i]]=MCMC(log(RW_init),p,logl_RW)
  
  #LGLL
  LGLL_opt = optim(par=c(0.0,log(RW_init)),fn=logl_LGLL,method="BFGS",control=list(fnscale=-1))
  res_LGLL[[i]]=MCMC(LGLL_opt$par,p,logl_LGLL)

}

##########################################################################################

saveRDS(res_RW, file="mcmc/monthly/MCMC_RW_20k.rds")
saveRDS(res_LGLL, file="mcmc/monthly/MCMC_LGLL_20k.rds")
