home=path.expand("~")
setwd(paste(home,"/CloudStation/UiS/Storage_model/",sep=""))
res_lin=readRDS(file="mcmc/monthly/MCMC_linear_20k.rds")
res_3=readRDS(file="mcmc/monthly/MCMC_RCS3_20k.rds")
res_7=readRDS(file="mcmc/monthly/MCMC_RCS7_20k.rds")
##################################################################
##################################################################
##################################################################
npar=10000
nburnin=2000
nMCMC = 10000+nburnin

dispnames=c("Natgas","Coffee","Cotton","Aluminum")
n=length(dispnames)
nparam=2

options(scipen=999)

delta_lin=matrix(NA,nrow=2,ncol=4)
delta_3=matrix(NA,nrow=2,ncol=4)
delta_7=matrix(NA,nrow=2,ncol=4)
b_lin=matrix(NA,nrow=2,ncol=4)
b_3=matrix(NA,nrow=2,ncol=4)
b_7=matrix(NA,nrow=2,ncol=4)

for (i in 1:n)
{
  t_lin=res_lin[[i]]$theta[(nburnin+1):nMCMC,3:4]
  t_lin[,1]=0.5+0.5*tanh(t_lin[,1])
  t_lin[,2]=exp(t_lin[,2])
  
  delta_lin[,i]=c(mean(t_lin[,1]),sd(t_lin[,1]))
  b_lin[,i]=c(mean(t_lin[,2]),sd(t_lin[,2]))
  #
  t_3=res_3[[i]]$theta[(nburnin+1):nMCMC,6:7]
  t_3[,1]=0.5+0.5*tanh(t_3[,1])
  t_3[,2]=exp(t_3[,2])
  
  delta_3[,i]=c(mean(t_3[,1]),sd(t_3[,1]))
  b_3[,i]=c(mean(t_3[,2]),sd(t_3[,2]))
  #
  t_7=res_7[[i]]$theta[(nburnin+1):nMCMC,10:11]
  t_7[,1]=0.5+0.5*tanh(t_7[,1])
  t_7[,2]=exp(t_7[,2])
  
  delta_7[,i]=c(mean(t_7[,1]),sd(t_7[,1]))
  b_7[,i]=c(mean(t_7[,2]),sd(t_7[,2]))
}
##########################################################################################################################
partxt=c("$\\delta$","$b$")
mtxt=c("Linear","RCS3","RCS7")
rowtxt=c(" ","Post. mean","Post. std.")

tab = matrix("",nrow=nparam*9,ncol=n+3)
tab[c(2,11),1]=partxt
tab[c(2,5,8,11,14,17),2]=rep(mtxt,2)
tab[,3]=rep(rowtxt,3*nparam)

specify_decimal = function(x, k) trimws(format(round(x, k), nsmall=k))

tab[2:3,4:7]=specify_decimal(delta_lin,4)
tab[5:6,4:7]=specify_decimal(delta_3,4)
tab[8:9,4:7]=specify_decimal(delta_7,4)

tab[c(11,14,17),4:7]=specify_decimal(rbind(b_lin[1,],b_3[1,],b_7[1,]),2)
tab[c(12,15,18),4:7]=specify_decimal(rbind(b_lin[2,],b_3[2,],b_7[2,]),3)

tab=rbind(c("","", "",dispnames),tab)
tab=tab[-2,]

library(xtable)


latextab=print(xtable(tab), hline.after = c(0,1,18),floating = FALSE,
               include.rownames = FALSE, include.colnames = FALSE,sanitize.text.function = function(x){x})

partab = paste("\\centering \n \\resizebox{0.65\\textwidth}{!}{ \n",latextab,"}\n",sep=" ")

cat(partab, file="Lyx/tab/partab_parametric.tex")
