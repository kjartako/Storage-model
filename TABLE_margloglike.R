C = 15

datnames=c("Natgas","Coffee","Cotton","Aluminum")
n=length(datnames)
##############
res_lin=readRDS(file=paste("mcmc/monthly/margloglike_linear_",C,".rds",sep=""))
res_RCS3=readRDS(file=paste("mcmc/monthly/margloglike_RCS3_",C,".rds",sep=""))
res_RCS7=readRDS(file=paste("mcmc/monthly/margloglike_RCS7_",C,".rds",sep=""))
res_LGLL=readRDS(file=paste("mcmc/monthly/margloglike_LGLL_",C,".rds",sep=""))

res_stoch=list()
for(i in 1:4)
{
  res_stoch[[i]]=readRDS(file=paste("mcmc/monthly/margloglike_stoch_",datnames[i],"_",C,".rds",sep=""))
}

##################################################################
##################################################################
##################################################################
options(scipen=999)
##########################################################################################################################
rtxt=c("Storage SSM","LGLL SSM"," ","Linear trend"," ","RCS3 trend"," ","RCS7 trend"," ")

marg_tab = matrix(NA,nrow=length(rtxt),ncol=n)

for(i in 1:n)
{
  marg_tab[1,i]=round(res_stoch[[i]]$marglike,2)
  
  marg_tab[2,i]=round(res_LGLL[[i]]$marglike,2)
  marg_tab[3,i]=round(marg_tab[1,i]-marg_tab[2,i],2)
  
  marg_tab[4,i]=round(res_lin[[i]]$marglike,2)
  marg_tab[5,i]=round(marg_tab[1,i]-marg_tab[4,i],2)
  
  marg_tab[6,i]=round(res_RCS3[[i]]$marglike,2)
  marg_tab[7,i]=round(marg_tab[1,i]-marg_tab[6,i],2)
  
  marg_tab[8,i]=round(res_RCS7[[i]]$marglike,2)
  marg_tab[9,i]=round(marg_tab[1,i]-marg_tab[8,i],2)
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

marg_tab=specify_decimal(marg_tab,2)
marg_tab[seq(3,9,2),] = paste0("(",marg_tab[seq(3,9,2),],")",sep="")

#####################################################
#####################################################
#####################################################

marg_tab=cbind(rtxt,marg_tab)
marg_tab=rbind(c(" ",datnames),marg_tab)

library(xtable)

latextab=print(xtable(marg_tab), hline.after = c(0,1,10),floating = FALSE,
               include.rownames = FALSE, include.colnames = FALSE)

partab = paste("\\centering \n \\resizebox{0.8\\textwidth}{!}{ \n",latextab,"}\n",sep=" ")

cat(partab, file="Lyx/tab/table_margloglike.tex")

########


test = readLines("Lyx/tab/table_margloglike.tex",-1)
test[9]=paste(test[9],"[0.3cm]",sep="")
test[11]=paste(test[11],"[0.2cm]",sep="")
test[13]=paste(test[13],"[0.2cm]",sep="")
test[15]=paste(test[15],"[0.2cm]",sep="")

writeLines(test,paste("Lyx/tab/table_margloglike_",C,".tex",sep=""))
