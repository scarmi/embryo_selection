
library(MASS)

risk_reduction_exclude = function(r2,K,q)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  integrand_exclude = function(t)
  {
    y = dnorm(t)*pnorm((zk-r*t)/sqrt(1-r^2),lower.tail=F)
    return(y)
  }
  risk = integrate(integrand_exclude,-Inf,zq)$value / (1-q)
  reduction = (K-risk)/K
  return(reduction)
}

risk_reduction_exclude2 = function(r2,K,q)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  a = zq*sqrt(2)
  integrand_t = function(t,x)
  {
    y = dnorm(t)*pnorm((zk-r/sqrt(2)*(x+t))/sqrt(1-r^2),lower.tail=F)
    return(y)
  }
  integrand_x = function(xs)
  {
    y = numeric(length(xs))
    for (i in seq_along(xs))
    {
      x = xs[i]
      temp = integrate(integrand_t,-Inf,a-x,x)$value
      denom = pnorm(a-x)
      if (denom==0) {denom=1e-300} # Avoid dividing by zero
      y[i] = dnorm(x)*temp/denom
    }
    return(y)
  }
  risk = integrate(integrand_x,-Inf,Inf)$value
  reduction = (K-risk)/K
  return(reduction)
}


r2 = 0.15
K = 0.01
n = 25
qs = seq(0,0.98,by=0.02)

nfam = 100000
scores = numeric(nfam)
disease_count = numeric(length(qs))
t_sick = qnorm(1-K)

for (i in seq_along(qs))
{
  cat('\r',i)
  q = qs[i]
  t_exclude = qnorm(1-q)*sqrt(r2)
  disease_count[i] = 0
  ws = rnorm(nfam,0,sqrt(r2/2))
  if (n<Inf)
  {
    xs = rnorm(nfam*n,0,sqrt(r2/2))
    xs = matrix(xs,nrow=nfam)
  }
  envs = rnorm(nfam,0,sqrt(1-r2))
  for (j in seq(1,nfam))
  {
    env = envs[j]
    w = ws[j]
    if (n==Inf)
    {
      while (T)
      {
       x = rnorm(1,0,sqrt(r2/2))
       score = x+w
       if (score<t_exclude)
         break
      }
      liab = score+env
    } else {
      x = xs[j,]
      scores = x+w
      if (any(scores<t_exclude))
      {
        s = min(which(scores<t_exclude))
      } else {
        s = 1
      }
      liab = scores[s]+env
    }
    if (liab>t_sick)
    {
      disease_count[i] = disease_count[i]+1
    }
  }
}
risk_red_sim = (K-disease_count/nfam)/K

risk_red_th = numeric(length(qs))
for (i in seq_along(qs))
{
  q = qs[i]
  risk_red_th[i] = risk_reduction_exclude2(r2,K,q)
}

#png("../Figures/ExcludeHigh.png",1400,1000,res=300)
par(mar=c(4.5,5,1,1))

plot(qs*100,risk_red_sim*100,cex=1.5,xlab="Percentile PS to exclude",ylab="Risk reduction(%)",cex.lab=1.5,cex.axis=1.25,lwd=2,main='Exclude high-risk',font.main=1,ylim=c(0,100))

lines(qs*100,risk_red_th*100,lwd=2,col='darkslategray',lty=1)
legend('bottomright',c(sprintf("Prevalence=%g",K),as.expression(bquote('r'^'2'*'='*.( r2)))),bty='n',cex=1.3,text.col='darkblue')

index = 2
segments(qs[index]*100,-10,qs[index]*100,risk_red_sim[index]*100,col='darkgray',lty=1,lwd=2)
segments(qs[index]*100,risk_red_sim[index]*100,-10,risk_red_sim[index]*100,col='darkgray',lty=1,lwd=2)

#dev.off()