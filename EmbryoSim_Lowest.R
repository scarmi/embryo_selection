
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

simulate_lowest_risk = function(r2,K,ns,nfam=10000)
{
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
  return(risk_red_sim)
}

