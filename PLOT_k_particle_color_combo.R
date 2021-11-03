home=path.expand("~")
setwd(paste(home,"/CloudStation/UiS/Storage_model/",sep=""))

library(RcppStorage)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
npar=10000
ndataset=4
a=0
C=10

datnames=c("Natural gas","Coffee","Cotton","Aluminum")
lineSize=0.75

kres = readRDS(file="mcmc/monthly/k_stoch.rds")

# k-plot
##############################################################################################

# Hexadecimal color specification 
Set1Col=brewer.pal(n = 3, name = "Set1")

brkpnts = matrix(c(37,97,157,217,rep(c(13,133,253,349),3)),nrow=4)
brklabs = matrix(c("2000","2005","2010","2015",rep(c("1990","2000","2010","2018"),3)),nrow=4)

Yscale = function(x) sprintf("%.1f", x)

plotlist=list()

for(i in 1:4)
{
  tmp=kres[[i]]
  startvals=tmp$startvals
  rep=tmp$rep
  N=length(tmp$p)
  
  kmean=apply(tmp$k,2,mean)
  k_low=apply(tmp$k, 2, quantile, probs = c(0.025))
  k_high=apply(tmp$k, 2, quantile, probs = c(0.975))
  xmean=apply(tmp$x,2,mean)
  x_low=apply(tmp$x, 2, quantile, probs = c(0.025))
  x_high=apply(tmp$x, 2, quantile, probs = c(0.975))
  
  pkmean=apply(tmp$pkdiff,2,mean)
  pk_low=apply(tmp$pkdiff, 2, quantile, probs = c(0.025))
  pk_high=apply(tmp$pkdiff, 2, quantile, probs = c(0.975))
  logf_left=log(evalf(rep$xleft,rep$xleft,rep$xright,rep$grid,rep$vals,a,exp(startvals[3]),C))
  logf_right=log(evalf(rep$xright,rep$xleft,rep$xright,rep$grid,rep$vals,a,exp(startvals[3]),C))
  
  dat=data.frame(x=1:length(kmean),k=kmean,lq=k_low,uq=k_high,p=tmp$p,xval=xmean,lx=x_low,ux=x_high,xleft=rep(rep$xleft,N),xright=rep(rep$xright,N),
                 pk=pkmean,lpk=pk_low,upk=pk_high,pkleft=rep(logf_left,N),pkright=rep(logf_right,N))
  
  p_plot=ggplot(dat, aes(x)) + 
    geom_ribbon(aes(ymin=lq, ymax=uq,fill="95 % CI, k"),alpha=0.3)+
    geom_line(aes(y = p,colour = "log p"),size=lineSize) + 
    geom_line(aes(y = k,colour = "k"),size=lineSize) + 
    scale_color_brewer(palette="Set1")+ggtitle(datnames[i])+
    scale_fill_manual(values=c("black"),breaks=c("95 % CI, k"))+
    theme_bw()+theme(legend.position="none",text = element_text(size=15),plot.title = element_text(size = 22))+
    ylab("")+xlab("t")+scale_x_continuous(breaks=brkpnts[,i],labels= brklabs[,i],expand=c(0.05, 0))+
    scale_y_continuous(labels=Yscale)
  
  pk_plot=ggplot(dat, aes(x)) + 
    geom_ribbon(aes(ymin=lpk, ymax=upk,fill="95 % CI, log f(x)"),alpha=0.3)+
    geom_line(aes(y = pk,colour = "log f(x)",linetype="log f(x)"),size=lineSize) + 
    geom_line(aes(y = pkleft,colour = "regime boundaries",linetype="regime boundaries"),size=lineSize)+
    geom_line(aes(y = pkright,colour = "regime boundaries",linetype="regime boundaries"),size=lineSize)+
    scale_colour_manual(values=c("log f(x)"=Set1Col[1],"regime boundaries"="gray20"),breaks=c("log f(x)","regime boundaries"))+
    scale_linetype_manual(values=c("solid","longdash"),breaks=c("log f(x)","regime boundaries"))+
    scale_fill_manual(values=c("black"),breaks=c("95 % CI, log f(x)"))+
    theme_bw()+theme(legend.position="none",text = element_text(size=15))+
    ylab("")+xlab("t")+scale_x_continuous(breaks=brkpnts[,i],labels= brklabs[,i],expand=c(0.05, 0))+
    scale_y_continuous(labels=Yscale)
  
  plotlist[[i]]=grid.arrange(p_plot,pk_plot)
}

emptyPlot = ggplot()+theme_void()

pdf(file="Lyx/fig/k_particle.pdf",width=2*7,height=3.5*4)
grid.arrange(plotlist[[1]],plotlist[[2]],emptyPlot,emptyPlot,plotlist[[3]],plotlist[[4]],nrow=3,heights=c(1,0.1,1))
dev.off()

##############################################################################################
##############################################################################################
##############################################################################################

# k-plot simulated data
##############################################################################################
simdata=readRDS(file="mcmc/monthly/simdata_natgas.rds")
ksim = readRDS(file="mcmc/monthly/k_stoch_sim.rds")

startvals=ksim$startvals
rep=ksim$rep
N=length(ksim$p)

kmean=apply(ksim$k,2,mean)
k_low=apply(ksim$k, 2, quantile, probs = c(0.025))
k_high=apply(ksim$k, 2, quantile, probs = c(0.975))
xmean=apply(ksim$x,2,mean)
x_low=apply(ksim$x, 2, quantile, probs = c(0.025))
x_high=apply(ksim$x, 2, quantile, probs = c(0.975))

pkmean=apply(ksim$pkdiff,2,mean)
pk_low=apply(ksim$pkdiff, 2, quantile, probs = c(0.025))
pk_high=apply(ksim$pkdiff, 2, quantile, probs = c(0.975))
logf_left=log(evalf(rep$xleft,rep$xleft,rep$xright,rep$grid,rep$vals,a,b=exp(startvals[3]),C))
logf_right=log(evalf(rep$xright,rep$xleft,rep$xright,rep$grid,rep$vals,a,b=exp(startvals[3]),C))

dat=data.frame(x=1:length(kmean),k=kmean,lq=k_low,uq=k_high,p=ksim$p,xval=xmean,lx=x_low,ux=x_high,xleft=rep(rep$xleft,N),xright=rep(rep$xright,N),
               pk=pkmean,lpk=pk_low,upk=pk_high,pkleft=rep(logf_left,N),pkright=rep(logf_right,N))

pdf(file="Lyx/fig/k_particle_sim.pdf",width=7,height=3.5*2)

p_plot=ggplot(dat, aes(x)) + 
  geom_ribbon(aes(ymin=lq, ymax=uq,fill="95 % CI, k"),alpha=0.3)+
  geom_line(aes(y = p,colour = "log p"),size=0.8*lineSize) + 
  geom_line(aes(y = k,colour = "k"),size=0.8*lineSize) + 
  geom_line(aes(y = simdata$k,colour = "real k"),size=0.8*lineSize) + 
  scale_color_brewer(palette="Set1")+
  scale_fill_manual(values=c("black"),breaks=c("95 % CI, k"))+
  theme_bw()+theme(legend.position="none",text = element_text(size=15))+ylab("")+xlab("t")

pk_plot=ggplot(dat, aes(x)) + 
  geom_ribbon(aes(ymin=lpk, ymax=upk,fill="95 % CI, log f(x)"),alpha=0.3)+
  geom_line(aes(y = pk,colour = "log f(x)",linetype="log f(x)"),size=lineSize) + 
  geom_line(aes(y = simdata$p-simdata$k,colour = "real log f(x)",linetype="real log f(x)"),size=lineSize) + 
  geom_line(aes(y = pkleft,colour = "regime boundaries",linetype="regime boundaries"),size=lineSize)+
  geom_line(aes(y = pkright,colour = "regime boundaries",linetype="regime boundaries"),size=lineSize)+
  scale_colour_manual(values=c("log f(x)"=Set1Col[1],"real log f(x)"=Set1Col[3],"regime boundaries"="gray20"),breaks=c("log f(x)","real log f(x)","regime boundaries"))+
  scale_linetype_manual(values=c("solid","solid","longdash"),breaks=c("log f(x)","real log f(x)","regime boundaries"))+
  scale_fill_manual(values=c("black"),breaks=c("95 % CI, log f(x)"))+
  theme_bw()+theme(legend.position="none",text = element_text(size=15))+ylab("")+xlab("t")

grid.arrange(p_plot,pk_plot)

dev.off()
