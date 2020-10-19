
library(MASS)

risk_reduction_lowest = function(r2,K,n)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  integrand_lowest = function(t)
  {
    arg = (zk-t*sqrt(1-r2/2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  risk = integrate(integrand_lowest,-Inf,Inf)$value
  reduction = (K-risk)/K
  return(reduction)
}

r2 = 0.15
K = 0.01
ns = seq(1,20)

nfam = 1000000
scores = numeric(nfam)
disease_count = numeric(length(ns))
t_sick = qnorm(1-K)

for (i in seq_along(ns))
{
  cat('\r',i)
  
  n = ns[i]
  disease_count[i] = 0
  ws = rnorm(nfam,0,sqrt(r2/2))
  xs = rnorm(nfam*n,0,sqrt(r2/2))
  xs = matrix(xs,nrow=nfam)
  envs = rnorm(nfam,0,sqrt(1-r2))
  for (j in seq(1,nfam))
  {
    w = ws[j]
    x = xs[j,]
    scores = x+w
    env = envs[j]
    liab = min(scores)+env
    if (liab>t_sick)
    {
      disease_count[i] = disease_count[i]+1
    }
  }
}
risk_red_sim = (K-disease_count/nfam)/K

risk_red_th = numeric(length(ns))
for (i in seq_along(ns))
{
  n = ns[i]
  risk_red_th[i] = risk_reduction_lowest(r2,K,n)
}

png("../Figures/LowestRisk.png",1400,1000,res=300)
par(mar=c(4.5,5,1,1))

plot(ns,risk_red_sim*100,cex=1.5,xlab="Number of embryos",ylab="Risk reduction(%)",cex.lab=1.5,cex.axis=1.25,lwd=2,main='Lowest risk',font.main=1,ylim=c(0,100))
lines(ns,risk_red_th*100,lwd=2,col='DarkSlateGray',lty=1)
legend('bottomright',c(sprintf("Prevalence=%g",K),as.expression(bquote('r'^'2'*'='*.( r2)))),bty='n',cex=1.3,text.col='darkblue')

index = 4
segments(ns[index],-10,ns[index],risk_red_sim[index]*100,col='darkgray',lty=2,lwd=2)
segments(ns[index],risk_red_sim[index]*100,-10,risk_red_sim[index]*100,col='darkgray',lty=2,lwd=2)

dev.off()