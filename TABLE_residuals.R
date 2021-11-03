tmp=readRDS(file="mcmc/monthly/RESIDUALS.rds")
PIT_residuals=tmp[[1]]
Pearson_residuals=tmp[[2]]

datnames=c("Natgas","Coffee","Cotton","Aluminum")
ndataset=length(datnames)
################################################################################################
Diagnostics=function(eps,eta,df)
{
  y1=skewness(eps)
  y2=kurtosis(eps)
  y3=jarque.bera.test(eps)$p.value
  y4=acf(eta,plot=FALSE)$acf[2]
  y5=Box.test(eta, lag = 12, type = c("Ljung-Box"), fitdf = df)$p.value
  y6=Box.test((eta^2), lag = 12, type = c("Ljung-Box"), fitdf = df)$p.value
  return(c(y1,y2,y3,y4,y5,y6))
}

################################################################################################
# Table generation
library(moments)
library(tseries)
library(fDMA)
nparams = c(3,4,7,11,1)
residual_names = " & Skew($\\xi_t$)   & Kurt($\\xi_t$)  & JB($\\xi_t$) & $\\rho_1(\\eta_t)$ & LB$_{12}(\\eta_t)$ & LB$_{12}(\\eta_t^2)$ \\\\"
################################################################################################
################################################################################################
stochastic = matrix(NA,ncol=6,nrow=ndataset)
linear = matrix(NA,ncol=6,nrow=ndataset)
RCS3 = matrix(NA,ncol=6,nrow=ndataset)
RCS7 = matrix(NA,ncol=6,nrow=ndataset)
LGLL = matrix(NA,ncol=6,nrow=ndataset)
for(j in 1:ndataset)
{
  stochastic[j,]=round(Diagnostics(eps=qnorm(PIT_residuals[[j]][1,]),eta=Pearson_residuals[[j]][1,],df=nparams[1]),3)
  linear[j,]=round(Diagnostics(eps=qnorm(PIT_residuals[[j]][2,]),eta=Pearson_residuals[[j]][2,],df=nparams[2]),3)
  RCS3[j,]=round(Diagnostics(eps=qnorm(PIT_residuals[[j]][3,]),eta=Pearson_residuals[[j]][3,],df=nparams[3]),3)
  RCS7[j,]=round(Diagnostics(eps=qnorm(PIT_residuals[[j]][4,]),eta=Pearson_residuals[[j]][4,],df=nparams[4]),3)
  LGLL[j,]=round(Diagnostics(eps=qnorm(PIT_residuals[[j]][5,]),eta=Pearson_residuals[[j]][5,],df=nparams[5]),3)
}

################################################################################################
################################################################################################
library(xtable)

stochtab=print(xtable(cbind(datnames,stochastic)),hline.after = c(), floating = FALSE, only.contents = TRUE,
             include.rownames = FALSE, include.colnames = FALSE,sanitize.text.function = function(x){x})
lintab=print(xtable(cbind(datnames,linear)),hline.after = c(),floating = FALSE, only.contents = TRUE,
            include.rownames = FALSE, include.colnames = FALSE,sanitize.text.function = function(x){x})
RCS3tab=print(xtable(cbind(datnames,RCS3)),hline.after = c(),floating = FALSE, only.contents = TRUE,
              include.rownames = FALSE, include.colnames = FALSE,sanitize.text.function = function(x){x})
RCS7tab=print(xtable(cbind(datnames,RCS7)),hline.after = c(),floating = FALSE, only.contents = TRUE,
              include.rownames = FALSE, include.colnames = FALSE,sanitize.text.function = function(x){x})
LGLLtab=print(xtable(cbind(datnames,LGLL)),hline.after = c(),floating = FALSE, only.contents = TRUE,
              include.rownames = FALSE, include.colnames = FALSE,sanitize.text.function = function(x){x})


################################################################################################

tab=paste("\\begin{tabular}{lcccccc}\n \\hline\\\\[-0.2cm] \n", residual_names,"\\hline\\\\[-0.2cm] \n",
          " & \\multicolumn{6}{c}{Storage SSM} \\\\ \n \\cline{2-7}\\\\ \n",stochtab,"[0.2cm]",
          " & \\multicolumn{6}{c}{LGLL SSM} \\\\ \n \\cline{2-7}\\\\ \n",LGLLtab,"[0.2cm]",
          " & \\multicolumn{6}{c}{Linear trend} \\\\ \n \\cline{2-7}\\\\ \n",lintab,"[0.2cm]",
          " & \\multicolumn{6}{c}{RCS3 trend} \\\\ \n \\cline{2-7}\\\\ \n",RCS3tab,"[0.2cm]",
          " & \\multicolumn{6}{c}{RCS7 trend} \\\\ \n \\cline{2-7}\\\\ \n",RCS7tab,"[0.2cm]",
          "\\hline \n \\end{tabular}")

tab = paste("\\centering\n\\resizebox{0.75\\textwidth}{!}{ \n",tab,"\n}",sep="")
tab = gsub("& 0 ","& $<$0.001 ",tab)
cat(tab, file="Lyx/tab/table_residuals.tex")


