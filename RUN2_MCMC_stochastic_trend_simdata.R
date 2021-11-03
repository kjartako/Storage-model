nburnin=2000
nMCMC = 10000+nburnin

natgas_MCMC=readRDS(file="mcmc/monthly/Natgas_10000.rds")$theta[(nburnin+1):nMCMC,]
postmean=colMeans(cbind(exp(natgas_MCMC[,1]),0.5+0.5*tanh(natgas_MCMC[,2]),exp(natgas_MCMC[,3])))
postmean_transformed=c(log(postmean[1]),atanh(2*postmean[2]-1),log(postmean[3]))

library(RcppStorage)
set.seed(1234)

a=0
C = 10
r = ((1.05)^(1/12))-1
v_x=1

ngrid = 20
nint = 128
int.grid = seq(from=-4.0,to=4.0,length.out=nint)
int.wts = dnorm(int.grid)*0.5*(int.grid[2]-int.grid[1])
xleft=0.0
xright=C
grid=seq(from=0,to=C,length.out=ngrid)
vals=seq(from=0,to=C,length.out=ngrid)
############################################################
# Particle marginal MH sampler
Nburnin=2000
Nburnin_adapt=500
Nmcmc = 10000+Nburnin
# number of particles
n = 10000
############################################################
simdata = readRDS("mcmc/monthly/simdata_natgas.rds")
p = simdata$p
N=length(p)

logl <- function(x){
  tryCatch(expr = SMC_sample(x,p,n,N,xleft,xright,grid,vals,a,C,r,v_x,int.grid,int.wts,nint,ngrid,30,0.5),
           ## But if an error occurs, do the following:
           error=function(error_message) {
             message(error_message)
             return(list(loglike=-5000,loglike_noprior=-5000))
           }
  )
}

PMCMC = function(startv)
{
  theta=matrix(rep(0,Nmcmc*3),ncol=3)
  #theta = c(log(v),delta,log(b))
  theta[1,] = startv
  
  loglike=numeric(Nmcmc)
  loglike_noprior=numeric(Nmcmc)
  acc = numeric(Nmcmc-1)
  accprobs = numeric(Nmcmc-1)
  
  tmp=logl(theta[1,])
  loglike[1]= tmp$loglike
  loglike_noprior[1]= tmp$loglike_noprior
  
  C_0 = diag(c(0.001,0.001,0.001))
  
  library(mvtnorm)
  
  for(t in 2:(Nburnin/2))
  {
    p_theta = rmvnorm(1,mean=theta[t-1,],sigma=C_0)
    
    loglt=logl(p_theta)
    p_loglike = loglt$loglike
    
    if(!is.finite(p_loglike))
    {
      p_loglike=-5000;
    }
    
    accp = exp(min(0,p_loglike-loglike[t-1]))
    
    accprobs[t-1]=accp
    
    if(runif(1)<accp)
    {
      theta[t,]=p_theta
      #xvals[t,]=p_x
      loglike[t]=p_loglike
      loglike_noprior[t]=loglt$loglike_noprior
      acc[t-1]=1
    }
    else
    {
      theta[t,]=theta[t-1,]
      #xvals[t,]=xvals[t-1,]
      loglike[t]=loglike[t-1]
      loglike_noprior[t]=loglike_noprior[t-1]
      acc[t-1]=0
    }
  }
  
  sd=(2.4^2)/3
  
  for(t in ((Nburnin/2)+1):Nburnin)
  {
    C_t = sd*cov(theta[Nburnin_adapt:(t-1),])+sd*0.000001*diag(3)
    p_theta = rmvnorm(1,mean=theta[t-1,],sigma=C_t)
    
    loglt=logl(p_theta)
    p_loglike = loglt$loglike
    
    if(!is.finite(p_loglike))
    {
      p_loglike=-5000;
    }
    
    accp = exp(min(0,p_loglike-loglike[t-1]))
    
    accprobs[t-1]=accp
    
    if(runif(1)<accp)
    {
      theta[t,]=p_theta
      #xvals[t,]=p_x
      loglike[t]=p_loglike
      loglike_noprior[t]=loglt$loglike_noprior
      acc[t-1]=1
    }
    else
    {
      theta[t,]=theta[t-1,]
      #xvals[t,]=xvals[t-1,]
      loglike[t]=loglike[t-1]
      loglike_noprior[t]=loglike_noprior[t-1]
      acc[t-1]=0
    }
  }
  
  for(t in (Nburnin+1):Nmcmc)
  {
    p_theta = rmvnorm(1,mean=theta[t-1,],sigma=C_t)
    
    loglt=logl(p_theta)
    p_loglike = loglt$loglike
    
    if(!is.finite(p_loglike))
    {
      p_loglike=-5000;
    }
    
    accp = exp(min(0,p_loglike-loglike[t-1]))
    
    accprobs[t-1]=accp
    
    if(runif(1)<accp)
    {
      theta[t,]=p_theta
      #xvals[t,]=p_x
      loglike[t]=p_loglike
      loglike_noprior[t]=loglt$loglike_noprior
      acc[t-1]=1
    }
    else
    {
      theta[t,]=theta[t-1,]
      #xvals[t,]=xvals[t-1,]
      loglike[t]=loglike[t-1]
      loglike_noprior[t]=loglike_noprior[t-1]
      acc[t-1]=0
    }
  }
  res = list(loglike=loglike,loglike_noprior=loglike_noprior,theta=theta,acc=acc,accprobs=accprobs,q_cov=C_t)
  
  return(res)
}

############################################################

res=PMCMC(postmean_transformed)

saveRDS(res, file=paste("mcmc/monthly/natgas_sim_",n,".rds",sep=""))
