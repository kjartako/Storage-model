library(readxl)
library(RcppStorage)

nburnin=2000
nMCMC = 10000+nburnin

dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])

datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])
datnames=c("Natgas","Coffee","Cotton","Aluminum") 

ndataset=length(datnames)

#########
a=0
C = 10
r = ((1.05)^(1/12))-1
v_x=1

##

ngrid = 20
nint = 128
int.grid = seq(from=-4.0,to=4.0,length.out=nint)
int.wts = dnorm(int.grid)*0.5*(int.grid[2]-int.grid[1])
xleft=0.0
xright=C
grid=seq(from=0,to=C,length.out=ngrid)
vals=seq(from=0,to=C,length.out=ngrid)

npar=10000
#######################################################################
#######################################################################
res=list()

for(i in 1:ndataset)
{
  p = log(datlist[[i]])
  N=length(p)
  
  theta=readRDS(file=paste("mcmc/monthly/",datnames[i],"_10000.rds",sep=""))$theta[(nburnin+1):nMCMC,]
  postmean=colMeans(cbind(exp(theta[,1]),0.5+0.5*tanh(theta[,2]),exp(theta[,3])))
  postmean_transformed=c(log(postmean[1]),atanh(2*postmean[2]-1),log(postmean[3]))
  
  set.seed(1234)
  tmp=SMC_sample(postmean_transformed,p,npar,N,xleft,xright,grid,vals,a,C,r,v_x,int.grid,int.wts,nint,ngrid,30,10)
  rep=tmp$rep
  x=tmp$Xp
  
  k=matrix(rep(NA,N*npar),nrow=npar)
  pkdiff=matrix(rep(NA,N*npar),nrow=npar)
  for(j in 1:N)
  {
    for (m in 1:npar)
    {
      k[m,j]=p[j]-log(evalf(x[m,j],rep$xleft,rep$xright,rep$grid,rep$vals,a,b=postmean[3],C))  
      pkdiff[m,j]=p[j]-k[m,j]
    }
  }
  res=c(res,list(list(name=datnames[i],p=p,rep=rep,x=x,k=k,pkdiff=pkdiff,startvals=postmean_transformed)))
}

saveRDS(res, file="mcmc/monthly/k_stoch.rds")
