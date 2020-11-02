
r2s = c(0.05,0.1,0.3)
Ks = c(0.01,0.05,0.2)
colors = c('Maroon','DarkCyan','DarkSlateBlue')

ns = seq(1,20)
labels = c('D','E','F')

png("../Figures/LowestRisk.png",2400,1000,res=300)
par(mfrow=c(1,length(Ks)))
par(mar=c(4,5.5,2,0.5))
risk_red_lowestrisk_sim_all = readRDS(file = './Data/select_lowest_risk_sim_maintext.rds')

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  plot(0,type="n",xlab="Number of embryos",ylab="Risk reduction (%)",cex.lab=1.5,cex.axis=1.25,main=sprintf('Lowest risk, K=%g',K),font.main=1,cex.main=1.5,ylim=c(0,100),xlim=c(1,20))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    risk_red_sim = risk_red_lowestrisk_sim_all[ki,ri,]
    risk_red_th = numeric(length(ns))
    for (i in seq_along(ns))
    {
      n = ns[i]
      risk_red_th[i] = risk_reduction_lowest(r2,K,n)
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
    # legs[ri] = as.expression(bquote('K='*.(K)*', r'^'2'*'='*.(r2)))
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }
  legend('bottomright',legs,col=colors,pch=0:2,bty='n',cex=1.3,text.col='DarkBlue',lwd=2)
  text(2,97,labels[ki],cex=2)
}
dev.off()


qs = seq(0,0.4,by=0.01)
n = 10
labels = c('A','B','C')
png("../Figures/ExcludeHigh_n10.png",2400,1000,res=300)
par(mfrow=c(1,length(Ks)))
par(mar=c(4,5.5,2,0.5))

risk_red_highrisk_sim_all = readRDS(file = "./Data/high_risk_exclude_sim_n10.rds")

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  plot(0,type='n',xlab="Percentile PS to exclude",ylab="Risk reduction(%)",cex.lab=1.5,cex.axis=1.25,cex.main=1.3,lwd=2,main=sprintf('High-risk exclusion, K=%g',K),font.main=1,ylim=c(0,100),xlim=c(0,max(qs)*100))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    risk_red_sim = risk_red_highrisk_sim_all[ki,ri,]
    risk_red_th = numeric(length(qs))
    for (i in seq_along(qs))
    {
      q = qs[i]
      risk_red_th[i] = risk_reduction_exclude(r2,K,q,n)
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
  text(4,97,labels[ki],cex=2)
}
dev.off()