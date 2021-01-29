
png('../Figures/DichotomizedTraits.png',2400,1200,res=300)

# IQ

r2 = 0.052
r = sqrt(r2)
K = pnorm(-2) # IQ<70, or 2SDs under the mean
ns = seq(1,20)
rrr = numeric(length(ns))

for (ni in seq_along(ns))
{
  n = ns[ni]
  rrr[ni] = risk_reduction_lowest(r2,K,n)
}

par(mar=c(4,5,2,1.5),mfrow=c(1,2))

plot(ns,rrr*100,lwd=2,xlab='Number of embryos',ylab='Relative risk reduction',cex.axis=1.5,cex.lab=1.5,main="IQ<70",font.main=1)

index = 5
segments(ns[index],-10,ns[index],rrr[index]*100,col='darkgrey',lty=1,lwd=2)
segments(ns[index],rrr[index]*100,-10,rrr[index]*100,col='darkgrey',lty=1,lwd=2)

text(15,10,'A',cex=1.5)

# LDL

r2 = 0.12
r = sqrt(r2)
Ks = seq(1,10)/100
n = 5
rrr = numeric(length(Ks))

for (Ki in seq_along(Ks))
{
  K = Ks[Ki]
  rrr[Ki] = risk_reduction_lowest(r2,K,n)
}

plot(Ks*100,rrr*100,lwd=2,xlab='LDL %-ile defined as affected',ylab='Relative risk reduction',cex.axis=1.5,cex.lab=1.25,main="LDL-C, n=5",font.main=1,ylim=c(0,70))

text(8,10,'B',cex=1.5)

dev.off()