
library(viridis)

r2 = 0.1
Ks = c(0.01,0.05,0.2)
rhos = c(-0.3,seq(-0.2,-0.05,by=0.05))

colors = magma(length(rhos)+1)[1:(length(rhos))]

input_file = "./Data/lowest_risk_two_taits.rds"
output_file = "../Figures/LowestRisk_TwoTraits.png"

ns = seq(1,20)
labels = c('A','B','C')
risk_red_lowest_sim_all = readRDS(file = input_file)
png(output_file,3000,1400,res=300)

par(mfrow=c(1,length(Ks)),mar=c(4,5.5,2,0.5),yaxs='i',xaxs='i')

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  plot(0,type='n',xlab="Number of embryos",ylab="Risk increase(%)",cex.lab=1.5,cex.axis=1.5,cex.main=1.3,lwd=2,,font.main=1,main=bquote(K==.(K)~","~r^2==.(r2)),ylim=c(0,100),xlim=c(0,max(ns)))
  
  #main=sprintf('K=%g',K)

  risk_red_target = risk_red_lowest_sim_all[ki,,1]
  points(ns,risk_red_target*100,cex=1.3,lwd=1.3,col='DarkSlateGray',pch=0,type='o')
  
  for (rhi in seq_along(rhos))
  {
    risk_red = risk_red_lowest_sim_all[ki,,rhi+1]
    points(ns,-risk_red*100,cex=1.3,lwd=1.3,col=colors[rhi],pch=1,type='o')
  }

  index = 5
  segments(ns[index],-10,ns[index],risk_red_target[index]*100,col='black',lty=1,lwd=1.5)
  segments(ns[index],risk_red_target[index]*100,-10,risk_red_target[index]*100,col='black',lty=1,lwd=1.5)
  
  legs = numeric()
  legs[1] = 'Selected trait (risk reduction)'
  for (rhi in seq_along(rhos))
  {
    legs[rhi+1] = as.expression(bquote("Correlated trait, "~rho==.(rhos[rhi])))
  }

  legend('topleft',legs,col=c('DarkSlateGray',colors),pch=c(0,rep(1,length(rhos))),lty=1,bty='n',cex=1.25)
  text(18,93,labels[ki],cex=1.8)
}
dev.off()
