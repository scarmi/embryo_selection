
prob_select = function(k,n,alpha)
{
  if (k==0)
    return(0)
  if (k==n)
    return(1)
  integrand = function(t)
  {
    y = k * dnorm(t,alpha,sqrt(1/2)) * pnorm(t,alpha,sqrt(1/2),lower.tail=F)^(k-1) * pnorm(t,0,sqrt(1/2),lower.tail=F)^(n-k)
    return(y)
  }
  return (integrate(integrand,-Inf,Inf)$value)
}

prob_unselect = function(k,n,alpha)
{
  integrand = function(t)
  {
    y = (n-k) * dnorm(t,0,sqrt(1/2)) * pnorm(t,0,sqrt(1/2),lower.tail=F)^(n-k-1) * pnorm(t,alpha,sqrt(1/2),lower.tail=F)^(k)
    return(y)
  }
  return (integrate(integrand,-Inf,Inf)$value)
}

new_freq = function(f,M,n,q)
{
  alpha = sqrt(1/(2*M*f*(1-f)))
  
  p = 0
  for (k in seq(0,n))
  {
    p = p + prob_select(k,n,alpha) * dbinom(k,n,0.5)
  }
  return (f^2 + 2*f*(1-f)*(q*p + (1-q)/2))
}

M = 500
n = 5
qs = c(1,2,5,10)/100

png("../Results/FreqChange.png",2000,1200,res=300)

library(viridis)
colors = viridis(length(qs))
# colors = c('darkblue','darkred','darkorchid','darkorange4')
plot(1,type="n",xlim=c(0,0.5),ylim=c(0,3),xlab='Initial MAF',ylab='MAF decrease (%)',cex.lab=1.5,cex.axis=1.5,main=sprintf('#embryos=%d, #loci=%d',n,M),font.main=1)
legends = numeric(length(qs))

for (qi in seq_along(qs))
{
  q = qs[qi]

  fs = seq(0.01,0.5,by=0.01)
  df = numeric(length(fs))
  for (fi in seq_along(fs))
  {
    f = fs[fi]
    fp = new_freq(f,M,n,q)
    df[fi] = 100*(f-fp)/f
  }
  
  points(fs,df,cex=1.25,col=colors[qi],lwd=2)
  legends[qi] = sprintf("Prop PES=%d%%",q*100)
}
legend("topright",legends,cex=1.25,col=colors,pch=1,pt.lwd=1.5,bty='n')

dev.off()