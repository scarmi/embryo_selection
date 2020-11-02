
r2s = c(0.05,0.1,0.3)
K = 0.05
qfs = c(50,98,98,98)
qms = c(50,25,50,98)

colors = c('Maroon','DarkCyan','DarkSlateBlue')

######### Exclude high risk ############

qs = seq(0,0.4,by=0.01)
n = 5
labels = c('A','B','C','D')
png("../Figures/ExcludeHigh_Conditional.png",3200,1000,res=300)
par(mfrow=c(1,length(qfs)))
par(mar=c(4,5.5,2,2))

risk_red_highrisk_sim_all = readRDS(file = "./Data/high_risk_exclude_sim_cond.rds")

for (qi in seq_along(qfs))
{
  qf = qfs[qi]
  qm = qms[qi]
  plot(0,type='n',xlab="Percentile PS to exclude",ylab="Risk reduction(%)",cex.lab=1.5,cex.axis=1.25,cex.main=1.15,lwd=2,main=sprintf('HR exclusion, par. 1=%g%%, par. 2=%g%%',qf,qm),font.main=1,ylim=c(0,100),xlim=c(0,max(qs)*100))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    risk_red_sim = risk_red_highrisk_sim_all[qi,ri,]
    risk_red_th = numeric(length(qs))
    for (i in seq_along(qs))
    {
      q = qs[i]
      risk_red_th[i] = risk_reduction_exclude_conditional(r2,K,q,n,qf,qm)
    }
    points(qs*100,risk_red_sim*100,cex=1.3,lwd=1.5,col=colors[ri],pch=ri-1)
    lines(qs*100,risk_red_th*100,lwd=1.8,col=colors[ri],lty=1)
    
    index = 3
    segments(qs[index]*100,-10,qs[index]*100,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
    segments(qs[index]*100,risk_red_sim[index]*100,-10,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
  }
  legs = numeric(length(r2s))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }
  pos = 'topright'
  legend(pos,legs,col=colors,pch=0:2,bty='n',cex=1.3,text.col='DarkBlue',lwd=2)
  text(4,97,labels[qi],cex=2)
}
dev.off()

######### Select lowest risk ############

ns = seq(1,20)
labels = c('E','F','G','H')

risk_red_lowestrisk_sim_all = readRDS(file = './Data/select_lowest_risk_sim_cond.rds')
png("../Figures/LowestRisk_Conditional.png",3200,1000,res=300)

par(mfrow=c(1,length(qfs)))
par(mar=c(4,5.5,2,2))

for (qi in seq_along(qfs))
{
  qf = qfs[qi]
  qm = qms[qi]
  plot(0,type="n",xlab="Number of embryos",ylab="Risk reduction (%)",cex.lab=1.5,cex.axis=1.25,main=sprintf('Lowest risk, par. 1=%g%%, par. 2=%d%%',qf,qm),font.main=1,cex.main=1.15,ylim=c(0,100),xlim=c(1,20))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    risk_red_sim = risk_red_lowestrisk_sim_all[qi,ri,]
    risk_red_th = numeric(length(ns))
    for (i in seq_along(ns))
    {
      n = ns[i]
      risk_red_th[i] = risk_reduction_lowest_conditional(r2,K,n,qf,qm)
    }
    points(ns,risk_red_sim*100,cex=1.3,lwd=1.5,col=colors[ri],pch=ri-1)
    lines(ns,risk_red_th*100,lwd=1.8,lty=1,col=colors[ri])
    
    index = 5
    segments(ns[index],-10,ns[index],risk_red_sim[index]*100,col='black',lty=1,lwd=2)
    segments(ns[index],risk_red_sim[index]*100,-10,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
  }
  legs = numeric(length(r2s))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }
  legend('bottomright',legs,col=colors,pch=0:2,bty='n',cex=1.3,text.col='DarkBlue',lwd=2)
  text(2,97,labels[qi],cex=2)
}
dev.off()
