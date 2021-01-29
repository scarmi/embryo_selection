
r2 = 0.1
r = sqrt(r2)

K = 0.05
zK = qnorm(1-K)

dr = 0.01
rs = seq(-2,2,by=dr)

pop_dist = dnorm(rs,0,r)
affected_dist = dnorm(rs,0,r) * pnorm((zK-rs)/sqrt(1-r2),lower.tail = F)/K

png('../Figures/ScoresDist.png',2000,1200,res=300)

par(mar=c(4,6,2,2))

plot(rs,affected_dist,type='l',lwd=2,col='darkred',xlab='PRS',ylab='Density',cex.axis=1.5,cex.lab=1.5)
lines(rs,pop_dist,lwd=2,col='blue')
legend('topright',c('Population','Affected'),bty = 'n',col=c('blue','darkred'),lwd=2)

dev.off()
