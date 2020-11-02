
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

risk_reduction_lowest_conditional = function(r2,K,n,qf,qm)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf/100)
  zqm = qnorm(qm/100)
  c = (zqf+zqm)/2 * r
  baseline= pnorm((zk-c)/sqrt(1-r^2/2),lower.tail=F)
  integrand_lowest_cond = function(t)
  {
    arg = (zk-c-t*sqrt(1-r^2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  risk = integrate(integrand_lowest_cond,-Inf,Inf)$value
  reduction = (baseline-risk)/baseline
  return(reduction)
}

simulate_lowest_risk = function(r2,K,ns,nfam=10000,parents_known=F,qf,qm)
{
  scores = numeric(nfam)
  disease_count = numeric(length(ns))
  disease_count_random = numeric(length(ns))
  baseline = numeric(length(ns))
  t_sick = qnorm(1-K)
  r = sqrt(r2)
  if (parents_known)
  {
    zqf = qnorm(qf/100)
    zqm = qnorm(qm/100)
    c = (zqf+zqm)/2 * r
  }
  for (i in seq_along(ns))
  {
    cat('\r',i)
    
    n = ns[i]
    disease_count[i] = 0
    cs = rnorm(nfam,0,r/sqrt(2))
    xs = rnorm(nfam*n,0,r/sqrt(2))
    xs = matrix(xs,nrow=nfam)
    envs = rnorm(nfam,0,sqrt(1-r2))
    for (j in seq(1,nfam))
    {
      if (!parents_known)
        c = cs[j]
      x = xs[j,]
      scores = x+c
      env = envs[j]
      liab = min(scores)+env
      if (liab>t_sick)
        disease_count[i] = disease_count[i]+1
      # Liability when not selecting
      liab = scores[1]+env
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