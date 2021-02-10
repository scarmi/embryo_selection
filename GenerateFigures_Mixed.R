
Ks = c(0.01,0.05,0.2)
colors = c('Maroon','DarkCyan','DarkSlateBlue')
r2 = 0.1

input_file_hre = "./Data/high_risk_exclude_sim_mixed.rds"
output_file_hre = "../Figures/ExcludeHigh_MixedStrategy.png"

qs = seq(0,1,by=0.02)
labels = c('A','B','C')
risk_red_highrisk_sim_all = readRDS(file = input_file_hre)
png(output_file_hre,2400,1000,res=300)

par(mfrow=c(1,length(Ks)),mar=c(4,5.5,2,0.5),yaxs='i',xaxs='i')

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  plot(0,type='n',xlab="Percentile PRS to exclude",ylab="Risk reduction(%)",cex.lab=1.5,cex.axis=1.25,cex.main=1.3,lwd=2,main=sprintf('K=%g',K),font.main=1,ylim=c(0,100),xlim=c(0,max(qs)*100))
  risk_red_sim_exclude = risk_red_highrisk_sim_all[ki,,1]
  risk_red_sim_mixed = risk_red_highrisk_sim_all[ki,,2]
  
  points(qs*100,risk_red_sim_exclude*100,cex=1.3,lwd=1.5,col=colors[1],pch=0)
  points(qs*100,risk_red_sim_mixed*100,cex=1.3,lwd=1.5,col=colors[2],pch=1)

  risk_red_th = numeric(length(qs))
  for (i in seq_along(qs))
  {
    q = qs[i]
    risk_red_th[i] = risk_reduction_exclude(r2,K,q,n)
  }
  lines(qs*100,risk_red_th*100,lwd=1.8,col=colors[1],lty=1)
  
  lowest_th = risk_reduction_lowest(r2,K,n)
  segments(0,lowest_th*100,100,lowest_th*100,col=colors[3],lty=2,lwd=1.75)
  
  #index = 3
  #segments(qs[index]*100,-10,qs[index]*100,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
  #segments(qs[index]*100,risk_red_sim[index]*100,-10,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
 
  legs = numeric()
  legs[1] = 'If all high-risk: choose random'
  legs[2] = 'If all high risk: choose lowest'
  legs[3] = 'Lowest risk prioritiztion'
  
  pos = 'topright'
  legend(pos,legs,col=colors,pch=c(0,1,NA),bty='n',cex=1,text.col='DarkBlue',lwd=c(2,NA,1.75),lty=c(1,1,2))
  text(9,75,labels[ki],cex=1.6)
}
dev.off()
