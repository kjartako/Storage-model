home=path.expand("~")
setwd(paste(home,"/CloudStation/UiS/Storage_model/",sep=""))

library(ggplot2)
library(gridExtra)
library(readxl)
library(RColorBrewer)

dat=as.matrix(read_excel("monthly_data.xlsx")[,c(2,9,12,15)])

datlist=list(natgas=na.omit(dat[,1]),coffee=dat[,2],cotton=dat[,3],alu=dat[,4])
datnames=c("Natural gas","Coffee","Cotton","Aluminum")
ndataset=length(datnames)

brkpnts = matrix(c(37,97,157,217,rep(c(13,133,253,349),3)),nrow=4)
brklabs = matrix(c("2000","2005","2010","2015",rep(c("1990","2000","2010","2018"),3)),nrow=4)

##################################################################################

k_stoch = readRDS(file="mcmc/monthly/k_stoch_smooth05.rds")
k_par = readRDS(file="mcmc/monthly/kvals_parametric.rds")

lineSize=0.8
ktypes=c("Stochastic","Linear","RCS3","RCS7")

plotlist=list()

for(i in 1:ndataset)
{
  p = log(datlist[[i]])
  N=length(p)
  kmat=k_par[[i]]
  
  plotdat=data.frame(x=1:N,y1=kmat[,1],y2=kmat[,2],y3=kmat[,3],y4=colMeans(k_stoch[[i]]$k))
  
  tmp_plot = ggplot(plotdat, aes(x)) + 
    geom_line(aes(y = y1,colour = "Linear",linetype="Linear",shape="Linear"),size=lineSize) + 
    geom_point(aes(y=y2,color = "RCS3",shape = "RCS3",linetype="RCS3"), size = 2,data=subset(plotdat, x %% 4 == 1))+
    geom_point(aes(y = y3,colour = "RCS7",linetype="RCS7",shape="RCS7"),size=2.5,data=subset(plotdat, x %% 4 == 1))+
    geom_line(aes(y = y4,colour = "Stochastic",linetype="Stochastic",shape="Stochastic"),size=lineSize) + 
    scale_linetype_manual(values=c("longdash","blank","blank","solid"),breaks=ktypes,name="test")+
    scale_shape_manual(values=c(32, 35, 0 ,32),breaks=ktypes,name="test")+
    theme_bw()+theme(legend.position="none",text = element_text(size=15),plot.title = element_text(size = 22))+
    ylab("")+xlab("t")+ggtitle(datnames[i])+guides(shape = guide_legend(override.aes = list(size = 2)))+
    scale_x_continuous(breaks=brkpnts[,i],labels= brklabs[,i],expand=c(0.05, 0))+scale_color_brewer(palette="Set1")
  tmp_plot
  
  plotlist[[i]]=tmp_plot
}

pdf("Lyx/fig/k_parametric.pdf",width=2*7,height=2*3.5)
parplot <- grid.arrange(arrangeGrob(plotlist[[1]] + theme(legend.position="none"),
                               plotlist[[2]] + theme(legend.position="none"),
                               plotlist[[3]] + theme(legend.position="none"),
                               plotlist[[4]] + theme(legend.position="none"),
                               nrow=2))
dev.off()


#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
# 
# mylegend<-g_legend(tmp_plot)
# 
# mylegend$grobs[[1]][[1]][[15]]$gp$lwd=2
# mylegend$grobs[[1]][[1]][[15]]$gp$fontsize=30
# mylegend$grobs[[1]][[1]][[16]]$gp$lwd=2
# mylegend$grobs[[1]][[1]][[16]]$gp$fontsize=30
# mylegend$grobs[[1]][[1]][[20]]$gp$lwd=2.5
# mylegend$grobs[[1]][[1]][[20]]$gp$fontsize=27
# mylegend$grobs[[1]][[1]][[21]]$gp$lwd=2.5
# mylegend$grobs[[1]][[1]][[21]]$gp$fontsize=27
# 
# 
# pdf("Lyx/fig/k_parametric.pdf",width=2*7,height=2*3.5)
# parplot <- grid.arrange(mylegend,
#                    arrangeGrob(plotlist[[1]] + theme(legend.position="none"),
#                                plotlist[[2]] + theme(legend.position="none"),
#                                plotlist[[3]] + theme(legend.position="none"),
#                                plotlist[[4]] + theme(legend.position="none"),
#                                nrow=2),
#                    nrow=2,heights=c(1,10))
# dev.off()

