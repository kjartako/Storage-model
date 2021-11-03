home=path.expand("~")
setwd(paste(home,"/CloudStation/UiS/Storage_model/",sep=""))
nburnin=2000
nMCMC = 10000+nburnin

natgas_MCMC=readRDS(file="mcmc/monthly/Natgas_10000.rds")$theta[(nburnin+1):nMCMC,]
postmean=colMeans(cbind(exp(natgas_MCMC[,1]),0.5+0.5*tanh(natgas_MCMC[,2]),exp(natgas_MCMC[,3])))

v=postmean[1]
delta = postmean[2]
b= postmean[3]
a=0

C = 10
r = ((1.05)^(1/12))-1
v_x = 1
beta = (1.0-delta)/(1.0+r)

ngrid = 20
nint = 128
int.grid = seq(from=-4.0,to=4.0,length.out=nint)
int.wts = dnorm(int.grid)*0.5*(int.grid[2]-int.grid[1])
############################################################
library(RcppStorage)
xleft=0.0
xright=C
grid=seq(from=0,to=C,length.out=ngrid)
vals=seq(from=0,to=C,length.out=ngrid)
rep=Storagefunc(xleft,xright,grid,vals,a,b,C,delta,beta,int.grid,int.wts,nint,ngrid,30)
############################################################
N_burnin = 500
N=300
N_tot = N_burnin+N

x = numeric(N)
k = numeric(N)
p = numeric(N)

set.seed(123)

x[1]=0
k[1]=0
p[1]=0

for(i in 2:N_tot)
{
  x[i]=((1-delta)*evalSigma(x[i-1],rep$xleft,rep$xright,rep$grid,rep$vals,C))+rnorm(1,0,v_x)
  k[i]=k[i-1]+rnorm(1,0,v)
  p[i]=k[i]+log(evalf(x[i],rep$xleft,rep$xright,rep$grid,rep$vals,a,b,C))
}

mydata = data.frame(x=x[(N_burnin+1):N_tot],k=k[(N_burnin+1):N_tot],p=p[(N_burnin+1):N_tot])

saveRDS(mydata,file="mcmc/monthly/simdata_natgas.rds")
