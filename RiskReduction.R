
library(MASS)

risk_reduction = function(r2,method,k,n=10,q=0.01)
{
  r = sqrt(r2)
  zk = qnorm(k, lower.tail=F)
  if (method==1) # Exclude high-risk, q-percentile
  {
    zq = qnorm(q, lower.tail=F)
    integrand1 = function(t)
    {
      y = dnorm(t)*pnorm((zk-r*t)/sqrt(1-r^2),lower.tail=F)
      return(y)
    }
    risk= integrate(integrand1,-Inf,zq)$value / (1-q)
  }
  if (method==2) # Select lowest-risk
  {
    integrand2 = function(t)
    {
      arg = (zk-t*sqrt(1-r2/2)) / (r/sqrt(2))
      y = dnorm(t)*pnorm(arg, lower.tail=F)^n
      return(y)
    }
    risk= integrate(integrand2,-Inf,Inf)$value
  }
  if (method==3) # Select high-risk (method 2)
  {
    integrand3 = function(t)
    {
      arg = (zk-t*sqrt(1-r2)) / r
      y = dnorm(t)*(1-pnorm(arg)/(1-q))
      return(y)
    }
    zq = qnorm(q, lower.tail=F)
    risk= integrate(integrand3,(zk-zq*r)/sqrt(1-r2),Inf)$value
  }
  if (method==4) # Select lowest-risk (method 2)
  {
    integrand4 = function(t)
    {
      arg = (zk-t*r/sqrt(2)) / sqrt(1-r2/2)
      y = n*dnorm(t)*pnorm(t,lower.tail = F)^(n-1)*pnorm(arg, lower.tail=F)
      return(y)
    }
    risk= integrate(integrand4,-Inf,Inf)$value
  }
  reduction = (k-risk)/k
  return(reduction)
}

r2 = 0.15
k = 0.01

qs = seq(0,0.99,by=0.01)
risk_red_1 = numeric(length(qs))
risk_red_3 = numeric(length(qs))

for (i in seq_along(qs))
{
  qi = qs[i]
  risk_red_1[i] = risk_reduction(r2,1,k,q=qi)
  risk_red_3[i] = risk_reduction(r2,3,k,q=qi)
}

ns = seq(1,20)
risk_red_2 = numeric(length(ns))
risk_red_4 = numeric(length(ns))

for (i in seq_along(ns))
{
  ni = ns[i]
  risk_red_2[i] = risk_reduction(r2,2,k,n=ni)
  risk_red_4[i] = risk_reduction(r2,4,k,n=ni)
}

# png("RiskReduction.png",2400,1500,res=300)
par(mfrow=c(1,2))

plot(qs*100,risk_red_1*100,type='l',cex=1.5,xlab="Percentile PS to exclude",ylab="Risk reduction(%)",cex.lab=1.5,cex.axis=1.25,lwd=2,main='Exclude high-risk',font.main=1)
index = 3
# segments(qs[index]*100,-10,qs[index]*100,risk_red_1[index]*100,col='darkgray',lty=1,lwd=2)
# segments(qs[index]*100,risk_red_1[index]*100,-10,risk_red_1[index]*100,col='darkgray',lty=1,lwd=2)

# plot(1/(1-qs),risk_red_1*100,type='l',cex=1.5,xlab="Effective number of embryos",ylab="Risk reduction(%)",cex.lab=1.5,cex.axis=1.25,lwd=2,main='Exclude high-risk',font.main=1,xlim=c(0,20))
legend('bottomright',c(sprintf("Prevalence=%g",k),as.expression(bquote('r'^'2'*'='*.( r2)))),bty='n',cex=1.3,text.col='darkblue')

plot(ns,risk_red_2*100,type='l',cex=1.5,xlab="Number of embryos",ylab="Risk reduction(%)",cex.lab=1.5,cex.axis=1.25,lwd=2,main='Select lowest-risk',font.main=1)

# dev.off()