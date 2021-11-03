C = 12.5


library(RcppStorage)
library(readxl)
library(splines)
dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])
datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])
##########################################################################################
##########################################################################################
##########################################################################################
n=4
linear=matrix(NA,nrow=n,ncol=4)
RCS3=matrix(NA,nrow=n,ncol=7)
RCS7=matrix(NA,nrow=n,ncol=11)

kpar = readRDS(file=paste("mcmc/monthly/k_par_optim_",C,".rds",sep=""))

for (i in 1:n)
{
  modpars = kpar[[i]]$pars
  # Linear
  linear[i,]=modpars[[1]]
  # RCS3
  RCS3[i,]=modpars[[6]]
  # RCS4
  RCS7[i,]=modpars[[11]]
}
##########################################################################################
##########################################################################################
##########################################################################################
logl <- function(kpars)
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
##########################################################################################
##########################################################################################
##########################################################################################
MCMC = function(kpars,p)
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
    C_t = sd*cov(theta[Nburnin_adapt:(t-1),])+sd*0.000001*diag(npar)
    
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
##########################################################################################
##########################################################################################
##########################################################################################
set.seed(1234)

a=0
r = ((1.05)^(1/12))-1
v_x=1

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
############################################################
# Particle marginal MH sampler
Nburnin=2000
Nburnin_adapt=500
Nmcmc = 20000+Nburnin
############################################################

res_lin=list()
res_3=list()
res_7=list()

for(j in 1:4)
{
  print(paste("data set",j,sep=" "))
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
  
  res_lin[[j]]=MCMC(linear[j,],p)
  #############################
  #############################
  #############################
  # RCS3
  knots=ns((1:N)/N, knots=quantile((1:N)/N,probs=c(0.25,0.5,0.75)),intercept = TRUE)
  
  k = function(t,kpars)
  {
    return(sum(knots[t,]*kpars))
  }
  
  res_3[[j]]=MCMC(RCS3[j,],p)
  #############################
  #############################
  #############################
  # RCS7
  knots=ns((1:N)/N, knots=quantile((1:N)/N,probs=c(0.125,0.250,0.375,0.500,0.625,0.750,0.875)),intercept = TRUE)
  
  k = function(t,kpars)
  {
    return(sum(knots[t,]*kpars))
  }
  res_7[[j]]=MCMC(RCS7[j,],p)
}


saveRDS(res_lin, file=paste("mcmc/monthly/MCMC_linear_20k_",C,".rds",sep=""))
saveRDS(res_3, file=paste("mcmc/monthly/MCMC_RCS3_20k_",C,".rds",sep=""))
saveRDS(res_7, file=paste("mcmc/monthly/MCMC_RCS7_20k_",C,".rds",sep=""))
