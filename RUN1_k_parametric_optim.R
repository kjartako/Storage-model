library(RcppStorage)
library(splines)
library(ggplot2)

C = 12.5

a=0
r = ((1.05)^(1/12))-1
v_x=1

c_alpha = 2
c_beta = 20

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
datnames=c("natgas","coffee","cotton","aluminium")
ndataset=length(datnames)

##

logl <- function(kpars)
{
  np=length(kpars)
  delta=0.5+0.5*tanh(kpars[np-1])
  b=exp(kpars[np])
  beta = (1.0-delta)/(1.0+r)
  
  rep= Storagefunc(xleft,xright,grid,vals,a,b, C, delta, beta,int.grid,int.wts,nint,ngrid, 30)
  xleft=rep$xleft
  xright=rep$xright
  grid=rep$grid
  vals=rep$vals
  
  x = numeric(N)
  x[1] = finv(exp(p[1]-k(1,kpars[1:(np-2)])),xleft,xright,grid,vals,a,b,C)
  
  loglike=0
  
  for(j in 2:N)
  {
    tmp=exp(p[j]-k(j,kpars[1:(np-2)]))
    x[j]=finv(tmp,xleft,xright,grid,vals,a,b,C)
    z =x[j]-(1-delta)*evalSigma(x[j-1],xleft,xright,grid,vals,C)
    Jacobian = abs(tmp/(evalf_dx(x[j], xleft, xright, grid, vals, a,b,C)))
    loglike = loglike + log((1/sqrt(2*pi))*exp(-0.5*z^2)*Jacobian)
  }
  
  delta_prior = log(delta)*(c_alpha-1.0)+log(1.0-delta)*(c_beta-1.0)+log(abs(-0.5+0.5*(2.0*delta-1.0)^2))
  b_prior = -0.5*(kpars[np])^2
  
  return(-loglike-delta_prior-b_prior)
}

##################################################################################
startvals=rbind(c(0.097,0.061,0.046,0.045),c(0.012,0.002,0.001,0.001),c(0.441,0.386,0.322,0.196))
startvals[1,]=log(startvals[1,])
startvals[2,]=atanh(2*startvals[2,]-1)
startvals[3,]=log(startvals[3,])
##################################################################################
##################################################################################
klist=list()

for(i in 1:ndataset)
{
  p = log(datlist[[i]])
  N=length(p)
  
  ###################################################################################################
  kmat=matrix(rep(0,N*3),ncol=3)
  
  ########################################
  ########################################
  ########################################
  # Linear
  
  k = function(t,kpars)
  {
    return(kpars[1]+kpars[2]*(t/N))
  }
  
  Par0=unname(c(lm(p ~ seq(from=0,to=1,length.out = N))$coefficients,startvals[2:3,i]))
  
  minval=logl(Par0)
  
  for(w in 1:5)
  {
    opt_lin=optim(Par0,logl)
    
    print(w)
    
    if(opt_lin$value<minval)
    {
      minval=opt_lin$value
      Par0=opt_lin$par
      print(opt_lin$par)
      print(opt_lin$value)
    }
  }
  
  for(q in 1:N)
  {
    kmat[q,1]=k(q,opt_lin$par)
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
  
  Par0=c(lm(p ~ knots-1)$coefficient,startvals[2:3,i])
  
  minval=logl(Par0)
  
  for(w in 1:5)
  {
    opt_3=optim(Par0,logl)
    
    print(w)
    
    if(opt_3$value<minval)
    {
      minval=opt_3$value
      Par0=opt_3$par
      print(opt_3$par)
      print(opt_3$value)
    }
  }
  
  for(q in 1:N)
  {
    kmat[q,2]=k(q,opt_3$par[1:5])
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
  
  Par0=c(lm(p ~ knots-1)$coefficient,startvals[2:3,i])
  
  minval=logl(Par0)
  
  for(w in 1:5)
  {
    opt_7=optim(Par0,logl)
    
    print(w)
    
    if(opt_7$value<minval)
    {
      minval=opt_7$value
      Par0=opt_7$par
      print(opt_7$par)
      print(opt_7$value)
    }
  }
  
  for(q in 1:N)
  {
    kmat[q,3]=k(q,opt_7$par[1:9])
  }
  ########################################
  ########################################
  ########################################
  ########################################
  
  klist[[i]]=list(kmat=kmat,pars=c(opt_lin,opt_3,opt_7))
}

saveRDS(klist,file=paste("mcmc/monthly/k_par_optim_",C,".rds",sep=""))


##################################################################################
##################################################################################
##################################################################################
