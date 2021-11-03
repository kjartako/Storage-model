C=12.5

npar=10000
nburnin=2000
nMCMC = 10000+nburnin

datnames=c("Natgas","Coffee","Cotton","Aluminum")
n=length(datnames)
nparam=3

library(LaplacesDemon)
options(scipen=999)

parmean=matrix(rep(0,n*3),nrow=n)
parsd=matrix(rep(0,n*3),nrow=n)
ESSvals=matrix(rep(0,n*3),nrow=n)

for (i in 1:n)
{
  res=readRDS(file=paste("mcmc/monthly/",datnames[i],"_",npar,"_",C,".rds",sep=""))$theta[(nburnin+1):nMCMC,]
  theta=cbind(exp(res[,1]),0.5+0.5*tanh(res[,2]),exp(res[,3])) 
  
  parmean[i,]=apply(theta,2,mean)
  parsd[i,]=apply(theta,2,sd)
  ESSvals[i,]=ESS(theta)
}
##########################################################################################################################
partxt=c("$v$","$\\delta$","$b$")
rowtxt=c(" ","Post. mean","Post. std.","ESS")

tab = matrix(rep(" ",(n+2)*nparam*4),nrow=nparam*4)

tab[,2]=rep(rowtxt,nparam)
tab[c(2,6,10),1]=partxt

for(i in 1:n)
{
  tab[c(2,6,10),i+2]=round(parmean[i,],4)
  tab[c(3,7,11),i+2]=round(parsd[i,],4)
  tab[c(4,8,12),i+2]=round(ESSvals[i,],0)
}

# tab[2,3:6]=round(as.numeric(tab[2,3:6]),3)
# tab[3,3:6]=round(as.numeric(tab[3,3:6]),4)
# tab[10,3:6]=round(as.numeric(tab[10,3:6]),2)
# tab[11,3:6]=round(as.numeric(tab[11,3:6]),3)

tab=rbind(c(" ", " ",datnames),tab)
tab=tab[-2,]

library(xtable)


latextab=print(xtable(tab), hline.after = c(0,1,12),floating = FALSE,
               include.rownames = FALSE, include.colnames = FALSE,sanitize.text.function = function(x){x})

tmptab = paste("\\centering \n \\resizebox{0.55\\textwidth}{!}{ \n",latextab,sep=" ")
partab = paste(tmptab,"}\n",sep=" ")

cat(partab, file=paste("Lyx/tab/partab_",C,".tex",sep=""))
