
nfam = 100000
r2 = 0.1
q = 0.4
t = qnorm(1-q)*sqrt(r2)
scores = numeric(nfam)

for (i in seq(1,nfam))
{
  parental_mean = rnorm(1,0,sqrt(r2/2))
  while(T)
  {
    embryo = rnorm(1,parental_mean,sqrt(r2/2))
    if (embryo<t)
    {
      scores[i] = embryo
      break
    }
  }
}

dx = 0.01
x = seq(-1.5,t+dx,by=dx)
hist(scores,breaks=x)
th = dnorm(x,0,sqrt(r2))/(1-q)*dx
lines(x,th*nfam,col='red',lwd=1,type='l')
