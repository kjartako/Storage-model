library(RcppStorage)

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

library(readxl)
dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])

datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])
datnames=c("Natgas","Coffee","Cotton","Aluminum")

ndataset=length(datnames)
npar=10000

fx=function(x,rep,b)
{
  return(evalf(x,rep$xleft,rep$xright,rep$grid,rep$vals,a,b,C))
}

fi=function(fval,rep,b)
{
  return(finv(fval, rep$xleft, rep$xright, rep$grid, rep$vals,a,b,C))
}

sig=function(x,rep)
{
  return(evalSigma(x, rep$xleft, rep$xright, rep$grid, rep$vals,C))
}

kalman.filter = function(sysq,sxsq,phi,mu1,omegasq,y)
{
  T = length(y)
  predmean = vector(length=T,mode="numeric")
  predvar = predmean
  filtmean = predmean
  filtvar = predmean
  loglikelihood = predmean
  resid = predmean

  PIT=numeric(T-1)
  Pearson=numeric(T-1)

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

    PIT[t-1]=pnorm(p[t],mean=predmean[t],sd=sqrt(predvar[t]+sysq))
    Pearson[t-1]=(p[t]-predmean[t])/sqrt(predvar[t]+sysq)
  }
  # return generated quantities
  return(list(predmean=predmean,predvar=predvar,filtmean=filtmean,
              filtvar=filtvar,loglikelihood=loglikelihood,resid=resid,PIT=PIT,Pearson=Pearson))
}
################################################################################################

nburnin=2000
nMCMC = 10000+nburnin

kpar=readRDS(file="mcmc/monthly/kvals_parametric.rds")

res_lin=readRDS(file="mcmc/monthly/MCMC_linear_20k.rds")
res_3=readRDS(file="mcmc/monthly/MCMC_RCS3_20k.rds")
res_7=readRDS(file="mcmc/monthly/MCMC_RCS7_20k.rds")
res_LGLL=readRDS(file="mcmc/monthly/MCMC_LGLL_20k.rds")
################################################################################################
################################################################################################
################################################################################################
PIT_residuals = list()
Pearson_residuals = list()

for(i in 1:ndataset)
{
  p = log(datlist[[i]])
  N=length(p)
  #################################################################
  # Stochastic
  set.seed(1234)

  stochMCMC=readRDS(file=paste("mcmc/monthly/",datnames[i],"_10000.rds",sep=""))$theta[(nburnin+1):nMCMC,]
  postmean=colMeans(cbind(exp(stochMCMC[,1]),0.5+0.5*tanh(stochMCMC[,2]),exp(stochMCMC[,3])))
  startvals=c(log(postmean[1]),atanh(2*postmean[2]-1),log(postmean[3]))

  filter=SMC_sample(startvals,p,npar,N,xleft,xright,grid,vals,a,C,r,v_x,int.grid,int.wts,nint,ngrid,30,0.5)
  rep=filter$rep
  x=filter$Xp

  v=postmean[1]
  delta=postmean[2]
  b=postmean[3]

  stochPearson=numeric(N-1)
  stochPIT=numeric(N-1)
  for(h in 2:N)
  {
    usum=0
    PearVec=numeric(npar)
    for(j in 1:npar)
    {
      xi = (1-delta)*sig(x[j,h-1],rep)+rnorm(1)
      usum=usum+pnorm(q=p[h],mean=p[h-1]+log(fx(xi,rep,b)/fx(x[j,h-1],rep,b)),sd=v)
      PearVec[j]=log(fx(xi,rep,b)/fx(x[j,h-1],rep,b))
    }
    stochPIT[h-1]=usum/npar
    stochPearson[h-1]=(p[h]-p[h-1]-mean(PearVec))/sqrt((v^2)+var(PearVec))
  }
  #################################################################
  # Linear
  k1=kpar[[i]][,1]

  t_lin=res_lin[[i]]$theta[(nburnin+1):nMCMC,3:4]
  t_lin=colMeans(cbind(0.5+0.5*tanh(t_lin[,1]),exp(t_lin[,2])))
  delta1=t_lin[1]
  b1=t_lin[2]
  beta1 = (1.0-delta1)/(1.0+r)
  linrep= Storagefunc(xleft,xright,grid,vals,a,b1, C, delta1, beta1 ,int.grid,int.wts,nint,ngrid, 30)

  linPearson=numeric(N-1)
  linPIT=numeric(N-1)
  for(h in 2:N)
  {
    linPIT[h-1]=pnorm(q=fi(exp(p[h]-k1[h]),linrep,b1)-(1-delta1)*sig(fi(exp(p[h-1]-k1[h-1]),linrep,b1),linrep),mean=0,sd=1)
    PearVec1=numeric(npar)
    for(j in 1:npar)
    {
      PearVec1[j]= k1[h]+ log(fx((1-delta1)*sig(fi(exp(p[h-1]-k1[h-1]),linrep,b1),linrep)+rnorm(1),linrep,b1))
    }
    linPearson[h-1]=(p[h]-mean(PearVec1))/sd(PearVec1)
  }
  #################################################################
  # RCS3
  k2=kpar[[i]][,2]

  t_3=res_3[[i]]$theta[(nburnin+1):nMCMC,6:7]
  t_3=colMeans(cbind(0.5+0.5*tanh(t_3[,1]),exp(t_3[,2])))
  delta2=t_3[1]
  b2=t_3[2]
  beta2 = (1.0-delta2)/(1.0+r)
  RCS3rep= Storagefunc(xleft,xright,grid,vals,a,b2, C, delta2, beta2 ,int.grid,int.wts,nint,ngrid, 30)

  RCS3Pearson=numeric(N-1)
  RCS3PIT=numeric(N-1)
  for(h in 2:N)
  {
    RCS3PIT[h-1]=pnorm(q=fi(exp(p[h]-k2[h]),RCS3rep,b2)-(1-delta2)*sig(fi(exp(p[h-1]-k2[h-1]),RCS3rep,b2),RCS3rep),mean=0,sd=1)
    PearVec2=numeric(npar)
    for(j in 1:npar)
    {
      PearVec2[j]= k2[h]+ log(fx((1-delta2)*sig(fi(exp(p[h-1]-k2[h-1]),RCS3rep,b2),RCS3rep)+rnorm(1),RCS3rep,b2))
    }
    RCS3Pearson[h-1]=(p[h]-mean(PearVec2))/sd(PearVec2)
  }
  #################################################################
  # RCS7
  k3=kpar[[i]][,3]

  t_7=res_7[[i]]$theta[(nburnin+1):nMCMC,10:11]
  t_7=colMeans(cbind(0.5+0.5*tanh(t_7[,1]),exp(t_7[,2])))
  delta3=t_7[1]
  b3=t_7[2]
  beta3 = (1.0-delta3)/(1.0+r)
  RCS7rep= Storagefunc(xleft,xright,grid,vals,a,b3, C, delta3, beta3 ,int.grid,int.wts,nint,ngrid, 30)

  RCS7Pearson=numeric(N-1)
  RCS7PIT=numeric(N-1)
  for(h in 2:N)
  {
    RCS7PIT[h-1]=pnorm(q=fi(exp(p[h]-k3[h]),RCS7rep,b3)-(1-delta3)*sig(fi(exp(p[h-1]-k3[h-1]),RCS7rep,b3),RCS7rep),mean=0,sd=1)
    PearVec3=numeric(npar)
    for(j in 1:npar)
    {
      PearVec3[j]= k3[h]+ log(fx((1-delta3)*sig(fi(exp(p[h-1]-k3[h-1]),RCS7rep,b3),RCS7rep)+rnorm(1),RCS7rep,b3))
    }
    RCS7Pearson[h-1]=(p[h]-mean(PearVec3))/sd(PearVec3)
  }
  #################################################################
  # LGLL
  t_LGLL = res_LGLL[[i]]$theta[(nburnin+1):nMCMC,]
  t_LGLL=colMeans(exp(t_LGLL))

  KalRes=kalman.filter(sysq=t_LGLL[1]^2,sxsq=t_LGLL[2]^2,phi=1.0,mu1=0.0,omegasq = 1000000.0,y=p)

  LGLL_Pearson=KalRes$Pearson
  LGLL_PIT=KalRes$PIT

  PIT_residuals[[i]] = rbind(stochPIT,linPIT,RCS3PIT,RCS7PIT,LGLL_PIT)
  Pearson_residuals[[i]] = rbind(stochPearson,linPearson,RCS3Pearson,RCS7Pearson,LGLL_Pearson)
}

saveRDS(list(PIT_residuals,Pearson_residuals),file="mcmc/monthly/RESIDUALS.rds")