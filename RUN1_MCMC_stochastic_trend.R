C = 9

set.seed(1234)
library(RcppStorage)

a=0
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
############################################################
# Particle marginal MH sampler
Nburnin=2000
Nburnin_adapt=500
Nmcmc = 10000+Nburnin
# number of particles
n = 10000
############################################################
############################################################
library(readxl)
dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])

datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])
datnames=c("Natgas","Coffee","Cotton","Aluminum")

startvals=rbind(c(0.097,0.061,0.046,0.045),c(0.012,0.002,0.001,0.001),c(0.441,0.386,0.322,0.196))
startvals[1,]=log(startvals[1,])
startvals[2,]=atanh(2*startvals[2,]-1)
startvals[3,]=log(startvals[3,])

############################################################
############################################################
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
      loglike[t]=p_loglike
      loglike_noprior[t]=loglt$loglike_noprior
      acc[t-1]=1
    }
    else
    {
      theta[t,]=theta[t-1,]
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
      loglike[t]=p_loglike
      loglike_noprior[t]=loglt$loglike_noprior
      acc[t-1]=1
    }
    else
    {
      theta[t,]=theta[t-1,]
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
      loglike[t]=p_loglike
      loglike_noprior[t]=loglt$loglike_noprior
      acc[t-1]=1
    }
    else
    {
      theta[t,]=theta[t-1,]
      loglike[t]=loglike[t-1]
      loglike_noprior[t]=loglike_noprior[t-1]
      acc[t-1]=0
    }
  }
  res = list(loglike=loglike,loglike_noprior=loglike_noprior,theta=theta,acc=acc,accprobs=accprobs,q_cov=C_t)
  
  return(res)
}

############################################################

for(i in 1:4)
{
  print(paste("data set",i,sep=" "))
  p = log(datlist[[i]])
  N=length(p)
  
  res=PMCMC(startvals[,i])
  
  saveRDS(res, file=paste("mcmc/monthly/",datnames[i],"_",n,"_",C,".rds",sep=""))  
}
