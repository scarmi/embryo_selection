
r2s = c(0.05,0.1,0.3)
Ks = c(0.01,0.05,0.2)
n = 5
colors = c('Maroon','DarkCyan','DarkSlateBlue')
rrr_breaks = seq(0,100,length.out=51)

input_file = "./Data/select_lowest_risk_personal_rrr_dist.rds"
output_file = "../Figures/LowestRisk_PersonalRRRDist.png"

labels = c('A','B','C')

risk_red_lowestrisk_th_all = readRDS(file = input_file)

#png(output_file,2400,1000,res=300)
par(mfrow=c(1,length(Ks)),mar=c(4,5.5,2,1),yaxs='i',xaxs='i')

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  plot(0,type="n",xlab="Relative risk reduction (%)",ylab="Probability",cex.lab=1.5,cex.axis=1.25,main=sprintf('K=%g',K,n),font.main=1,cex.main=1.5,ylim=c(0,0.5),xlim=c(0,100))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    risk_red = risk_red_lowestrisk_th_all[ki,ri,]
    rrr_dist = hist(risk_red*100,breaks=rrr_breaks,plot=F)
    rrr_dist_count = rrr_dist$counts/length(risk_red)
    rrr_dist_mids = rrr_dist$mids
    
    lines(rrr_dist_mids,rrr_dist_count,lwd=2,col=colors[ri])
    
    mean_rrr = mean(risk_red)
    mean_th_rrr = risk_reduction_lowest(r2,K,n)
    cat(sprintf("K=%g, r2=%g;  \tmean personal RRR=%.2f mean population RRR=%.2f\n",K,r2,mean_rrr,mean_th_rrr))
  }
  legs = numeric(length(r2s))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }
  legend('topright',legs,col=colors,pch=0:2,bty='n',cex=1.3,text.col='DarkBlue',lwd=2)
  text(9,0.45,labels[ki],cex=2)
}
dev.off()

output_file = "../Figures/LowestRisk_RRRvs_c.png"

png(output_file,2400,1000,res=300)
par(mfrow=c(1,length(Ks)),mar=c(4,5.5,2,1),yaxs='i',xaxs='i')
labels = c('D','E','F')
num_cs_bins = 10000

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  plot(0,type="n",xlab="Av parental PRS quantile",ylab="Relative risk reduction (%)",cex.lab=1.5,cex.axis=1.25,main=sprintf('K=%g',K,n),font.main=1,cex.main=1.5,ylim=c(0,100),xlim=c(0,100))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    risk_red = risk_red_lowestrisk_th_all[ki,ri,]
    
    qs = seq(0,100,length.out = num_cs_bins+1)
    dq = qs[2]-qs[1]
    qs = qs[1:(length(qs)-1)]+dq/2

    lines(qs,risk_red*100,lwd=2,col=colors[ri])
  }
  legs = numeric(length(r2s))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }
  legend('bottomleft',legs,col=colors,pch=0:2,bty='n',cex=1.3,text.col='DarkBlue',lwd=2)
  text(88,93,labels[ki],cex=2)
}
dev.off()