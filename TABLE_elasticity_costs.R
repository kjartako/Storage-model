C=12.5

set.seed(123)
library(readxl)
library(RcppStorage)

a=0
r = ((1.05)^(1/12))-1

ngrid = 20
nint = 128
int.grid = seq(from=-4.0,to=4.0,length.out=nint)
int.wts = dnorm(int.grid)*0.5*(int.grid[2]-int.grid[1])
xleft=0.0
xright=C
grid=seq(from=0,to=C,length.out=ngrid)
vals=seq(from=0,to=C,length.out=ngrid)

dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])

datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])
datnames=c("Natgas","Coffee","Cotton","Aluminum")

ndataset=length(datnames)


#res_lin=readRDS(file="mcmc/monthly/MCMC_linear_20k.rds")
res_3=readRDS(file=paste("mcmc/monthly/MCMC_RCS3_20k_",C,".rds",sep=""))
res_7=readRDS(file=paste("mcmc/monthly/MCMC_RCS7_20k_",C,".rds",sep=""))

kstoch = readRDS(file=paste("mcmc/monthly/k_stoch_smooth05_",C,".rds",sep=""))
kpar = readRDS(file=paste("mcmc/monthly/kvals_parametric_",C,".rds",sep=""))

nburnin=2000
nMCMC = 10000+nburnin
##############################################
##############################################
elast=matrix(rep(0,ndataset*3),ncol=ndataset)
cost=matrix(rep(0,ndataset*3),ncol=ndataset)

for(i in 1:ndataset)
{
  p = log(datlist[[i]])
  N=length(p)
  #
  kvals = kpar[[i]]
  #
  #t_lin=res_lin[[i]]$theta[(nburnin+1):nMCMC,3:4]
  #t_lin=colMeans(cbind(0.5+0.5*tanh(t_lin[,1]),exp(t_lin[,2])))
  #
  t_3=res_3[[i]]$theta[(nburnin+1):nMCMC,6:7]
  t_3=colMeans(cbind(0.5+0.5*tanh(t_3[,1]),exp(t_3[,2])))
  #
  t_7=res_7[[i]]$theta[(nburnin+1):nMCMC,10:11]
  t_7=colMeans(cbind(0.5+0.5*tanh(t_7[,1]),exp(t_7[,2])))
  #
  deltapar=c(t_3[1],t_7[1])
  bpar=c(t_3[2],t_7[2])

  # Stochastic trend
  stochMCMC=readRDS(file=paste("mcmc/monthly/",datnames[i],"_10000_",C,".rds",sep=""))$theta[(nburnin+1):nMCMC,]
  postmean=colMeans(cbind(exp(stochMCMC[,1]),0.5+0.5*tanh(stochMCMC[,2]),exp(stochMCMC[,3])))

  xvals = kstoch[[i]]$x
  b=postmean[3]
  delta=postmean[2]
  xmean = mean(xvals)
  elast[1,i]= -1.0/(b*xmean)
  cost[1,i]=-((1-delta)^12-1)

  # Parametric trends
  for(j in 1:2)
  {
    delta=deltapar[j]
    b = bpar[j]
    beta = (1.0-delta)/(1.0+r)

    rep= Storagefunc(xleft,xright,grid,vals,a,b, C, delta, beta,int.grid,int.wts,nint,ngrid, 30)

    xvals=numeric(N)
    for(q in 1:N)
    {
      xvals[q] = finv(exp(p[q]-kvals[q,j]),rep$xleft,rep$xright,rep$grid,rep$vals,a,b,C)
    }
    xmean = mean(xvals)
    elast[j+1,i]= -1.0/(b*xmean)
    cost[j+1,i]= -((1-delta)^12-1)
  }
}


##############################################
##############################################
##############################################
##############################################


specify_decimal = function(x, k) trimws(format(round(x, k), nsmall=k))
cost=specify_decimal(100*cost,1)
elast=specify_decimal(elast,2)

rnames = c("Storage SSM","RCS3 trend","RCS7 trend")


tabdat = cbind(rnames,cost[,1],elast[,1],rep(" ",3),cost[,2],elast[,2],rep(" ",3),cost[,3],elast[,3],rep(" ",3),cost[,4],elast[,4])

library(xtable)


tab_elascost=print(xtable(tabdat), hline.after = c(),floating = FALSE, only.contents = TRUE,
               include.rownames = FALSE, include.colnames = FALSE,sanitize.text.function = function(x){x})


elastab=paste("\\begin{tabular}{lrrcrrcrrcrr}\n \\hline \n",
          "&\\multicolumn{2}{c}{Natgas}  && \\multicolumn{2}{c}{Coffee} && \\multicolumn{2}{c}{Cotton} && \\multicolumn{2}{c}{Aluminum}\\\\ \n",
          "\\cline{2-3} \\cline{5-6} \\cline{8-9} \\cline{11-12} \\\\[-0.2cm] \n",
          "&\\multicolumn{1}{c}{costs} & \\multicolumn{1}{c}{elast.} && \\multicolumn{1}{c}{costs} & \\multicolumn{1}{c}{elast.} && \\multicolumn{1}{c}{costs} & \\multicolumn{1}{c}{elast.} && \\multicolumn{1}{c}{costs} & \\multicolumn{1}{c}{elast.}\\\\ \n",
          "\\hline\\\\[-0.2cm] \n",tab_elascost,
          "\\hline \n \\end{tabular}")


elastab = paste("\\centering \n \\resizebox{0.85\\textwidth}{!}{ \n",elastab,"}\n",sep=" ")

cat(elastab, file=paste("Lyx/tab/table_ela_costs_",C,".tex",sep=""))
  