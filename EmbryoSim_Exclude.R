
library(MASS)

risk_reduction_exclude_wrong = function(r2,K,q)
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

risk_reduction_exclude = function(r2,K,q)
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
      internal_int = integrate(integrand_t,-Inf,a-x,x)$value
      denom = pnorm(a-x)
      if (denom==0) {denom=1e-300} # Avoid dividing by zero
      y[i] = dnorm(x)*internal_int/denom
    }
    return(y)
  }
  risk = integrate(integrand_x,-Inf,Inf)$value
  reduction = (K-risk)/K
  return(reduction)
}

simulate_exclude_high = function(r2,K,n,qs,nfam=10000)
{
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
  return(risk_red_sim)
}