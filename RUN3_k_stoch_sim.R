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

npar=10000
nrep=100
#######################################################################
#######################################################################

theta=readRDS(file="mcmc/monthly/natgas_sim_10000.rds")$theta[2001:12000,]
parvals=c(mean(exp(theta[,1])),mean(0.5+0.5*tanh(theta[,2])),mean(exp(theta[,3])))
startvals=c(log(parvals[1]),atanh(2*parvals[2]-1),log(parvals[3]))


simdata = readRDS("mcmc/monthly/simdata_natgas.rds")
p = simdata$p
N=length(p)

set.seed(123)
tmp=SMC_sample(startvals,p,npar,N,xleft,xright,grid,vals,a,C,r,v_x,int.grid,int.wts,nint,ngrid,30,10)
rep=tmp$rep
x=tmp$Xp

k=matrix(rep(NA,N*npar),nrow=npar)
pkdiff=matrix(rep(NA,N*npar),nrow=npar)
for(j in 1:N)
{
  for (m in 1:npar)
  {
    k[m,j]=p[j]-log(evalf(x[m,j],rep$xleft,rep$xright,rep$grid,rep$vals,a,b=parvals[3],C))  
    pkdiff[m,j]=p[j]-k[m,j]
  }
}
res=list(name="sim",p=p,rep=rep,x=x,k=k,pkdiff=pkdiff,startvals=startvals)


saveRDS(res, file="mcmc/monthly/k_stoch_sim.rds")
