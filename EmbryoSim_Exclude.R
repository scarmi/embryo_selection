
library(MASS)

risk_reduction_exclude_appx1 = function(r2,K,q)
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

risk_reduction_exclude_appx2 = function(r2,K,q)
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

risk_reduction_exclude = function(r2,K,q,n)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  integrand_t = function(t,u)
  {
    y = dnorm(t)*pnorm((zk-r/sqrt(2)*(u+t))/sqrt(1-r^2),lower.tail=F)
    return(y)
  }
  integrand_u = function(us)
  {
    y = numeric(length(us))
    for (i in seq_along(us))
    {
      u = us[i]
      beta = zq*sqrt(2)-u
      internal_int = integrate(integrand_t,-Inf,beta,u)$value
      denom = pnorm(beta)
      if (denom==0) {denom=1e-300} # Avoid dividing by zero
      numer = dnorm(u)*(1-pnorm(beta,lower.tail=F)^n) * internal_int
      term1 = numer/denom
      
      internal_int = integrate(integrand_t,beta,Inf,u)$value
      denom = pnorm(beta,lower.tail=F)
      if (denom==0) {denom=1e-300} # Avoid dividing by zero
      numer = dnorm(u)*(pnorm(beta,lower.tail=F)^n) * internal_int
      term2 = numer/denom
      y[i] = term1 + term2
    }
    return(y)
  }
  risk = integrate(integrand_u,-Inf,Inf)$value
  reduction = (K-risk)/K
  return(reduction)
}

risk_reduction_exclude_conditional = function(r2,K,q,n,qf,qm)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  zqf = qnorm(qf/100)
  zqm = qnorm(qm/100)
  c = (zqf+zqm)/2 * r
  baseline= pnorm((zk-c)/sqrt(1-r^2/2),lower.tail=F)
  gamma = zq*sqrt(2) - c/(r/sqrt(2))
  
  integrand_t = function(t)
  {
    y = dnorm(t)*pnorm((zk-t*r/sqrt(2)-c)/sqrt(1-r^2),lower.tail=F)
    return(y)
  }

  internal_int = integrate(integrand_t,-Inf,gamma)$value
  denom = pnorm(gamma)
  numer = (1-pnorm(gamma,lower.tail=F)^n) * internal_int
  term1 = numer/denom
  
  internal_int = integrate(integrand_t,gamma,Inf)$value
  denom = pnorm(gamma,lower.tail=F)
  numer = pnorm(gamma,lower.tail=F)^n * internal_int
  term2 = numer/denom
    
  risk = term1 + term2
  
  reduction = (baseline-risk)/baseline
  return(reduction)
}

simulate_exclude_high = function(r2,K,n,qs,nfam=10000,parents_known,qf,qm)
{
  r = sqrt(r2)
  scores = numeric(nfam)
  disease_count = numeric(length(qs))
  disease_count_random = numeric(length(qs))
  baseline = numeric(length(qs))
  t_sick = qnorm(1-K)
  if (parents_known)
  {
    zqf = qnorm(qf/100)
    zqm = qnorm(qm/100)
    c = (zqf+zqm)/2 * r
  }
  for (i in seq_along(qs))
  {
    cat('\r',i)
    q = qs[i]
    t_exclude = qnorm(1-q)*sqrt(r2)
    disease_count[i] = 0
    cs = rnorm(nfam,0,sqrt(r2/2))
    if (n<Inf)
    {
      xs = rnorm(nfam*n,0,sqrt(r2/2))
      xs = matrix(xs,nrow=nfam)
    }
    envs = rnorm(nfam,0,sqrt(1-r2))
    for (j in seq(1,nfam))
    {
      env = envs[j]
      if (!parents_known)
        c = cs[j]
      if (n==Inf)
      {
        while (T)
        {
         x = rnorm(1,0,sqrt(r2/2))
         score = x+c
         if (score<t_exclude)
           break
        }
        liab = score+env
      } else {
        x = xs[j,]
        scores = x+c
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
      # Liability when not selecting
      x = rnorm(1,0,sqrt(r2/2))
      liab = x + c + env
      if (liab>t_sick)
        disease_count_random[i] = disease_count_random[i]+1
    }
    
    if (parents_known) {
      baseline[i] = disease_count_random[i]/nfam
    } else {
      baseline[i] = K
    }
  }
  risk_red_sim = (baseline-disease_count/nfam)/baseline
  return(risk_red_sim)
}