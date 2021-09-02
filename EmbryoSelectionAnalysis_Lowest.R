
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

risk_reduction_lowest_conditional = function(r2,K,n,qf,qm,relative=T,parental_avg_given=F)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf/100)
  zqm = qnorm(qm/100)
  if (parental_avg_given)
  {
    # It is assumed that what is given is directly the parental average, so that the paternal and maternal quantiles are the same (both representing the quantile of the parental average)
    c = zqf * r/sqrt(2)
  } else {
    c = (zqf+zqm)/2 * r
  }
  baseline = pnorm((zk-c)/sqrt(1-r^2/2),lower.tail=F)
  integrand_lowest_cond = function(t)
  {
    arg = (zk-c-t*sqrt(1-r^2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  risk = integrate(integrand_lowest_cond,-Inf,Inf,rel.tol=1e-9)$value
  if (relative) {
    reduction = (baseline-risk)/baseline
  } else {
    reduction = baseline-risk
  }
  return(reduction)
}

risk_reduction_lowest_family_history = function(r2,h2,K,n,df,dm)
{
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail=F)
  
  integrand_lowest_given_parents = function(t,gm,gf)
  {
    arg = (zk - (gm+gf)/2 - t*r/sqrt(2)) / sqrt(1-h^2/2-r^2/2)
    y = n * dnorm(t)*pnorm(t,lower.tail=F)^(n-1) * pnorm(arg, lower.tail=F)
    return(y)
  }
  
  posterior = function(gm,gf)
  {
    y = 1
    y = y * dnorm(gm/h)/h
    y = y * dnorm(gf/h)/h
    arg = (zk-gm)/sqrt(1-h^2)
    if (dm)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    arg = (zk-gf)/sqrt(1-h^2)
    if (df)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    return(y)
  }
  
  integrand_gm = function(gms,gf,baseline)
  {
    y = numeric(length(gms))
    for (i in seq_along(gms))
    {
      gm = gms[i]
      if (baseline)
      {
        arg = (zk - (gm+gf)/2) / sqrt(1-h^2/2)
        y[i] = pnorm(arg, lower.tail=F)
      } else {
        y[i] = integrate(integrand_lowest_given_parents,-Inf,Inf,gm,gf)$value
      }
      post = posterior(gm,gf)
      y[i] = y[i] * post
    }
    return(y)
  }
  
  integrand_gf = function(gfs,baseline)
  {
    y = numeric(length(gfs))
    for (i in seq_along(gfs))
    {
      gf = gfs[i]
      y[i] = integrate(integrand_gm,-Inf,Inf,gf,baseline)$value
    }
    return(y)
  }
  start_time = Sys.time()
  
  risk_selection = integrate(integrand_gf,-Inf,Inf,F)$value
  risk_baseline = integrate(integrand_gf,-Inf,Inf,T)$value
  
  end_time = Sys.time()
  cat(sprintf("Total time: %g\n",end_time-start_time))
  
  relative_reduction = (risk_baseline-risk_selection)/risk_baseline
  abs_reduction = risk_baseline-risk_selection
  return(c(risk_baseline,risk_selection,relative_reduction,abs_reduction))
}
  
# A slower implementation - do not use
risk_reduction_lowest_family_history_long = function(r2,h2,K,n,df,dm,relative=T)
{
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail=F)
  
  integrand_lowest_given_parents = function(t,sm,sf,gm,gf)
  {
    arg = (zk - (sm+gm+sf+gf)/2 - t*r/sqrt(2)) / sqrt(1-h^2/2-r^2/2)
    y = n * dnorm(t)*pnorm(t,lower.tail=F)^(n-1) * pnorm(arg, lower.tail=F)
    return(y)
  }
 
  posterior_all = function(sm,sf,gm,gf)
  {
    y = 1
    y = y * dnorm(sm/r)/r * dnorm(gm/sqrt(h^2-r^2))/sqrt(h^2-r^2)
    y = y * dnorm(sf/r)/r * dnorm(gf/sqrt(h^2-r^2))/sqrt(h^2-r^2)
    arg = (zk-sm-gm)/sqrt(1-h^2)
    if (dm)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    arg = (zk-sf-gf)/sqrt(1-h^2)
    if (df)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    return(y)
  }
   
  integrand_sm = function(sms,sf,gm,gf,baseline)
  {
    y = numeric(length(sms))
    for (i in seq_along(sms))
    {
      sm = sms[i]
      if (baseline)
      {
        arg = (zk - (sm+gm+sf+gf)/2) / sqrt(1-h^2/2)
        y[i] = pnorm(arg, lower.tail=F)
      } else {
        y[i] = integrate(integrand_lowest_given_parents,-Inf,Inf,sm,sf,gm,gf)$value
      }
      post = posterior_all(sm,sf,gm,gf)
      y[i] = y[i] * post
    }
    return(y)
  }
  integrand_sf = function(sfs,gm,gf,baseline)
  {
    y = numeric(length(sfs))
    for (i in seq_along(sfs))
    {
      sf = sfs[i]
      y[i] = integrate(integrand_sm,-Inf,Inf,sf,gm,gf,baseline)$value
    }
    return(y)
  }
  integrand_gm = function(gms,gf,baseline)
  {
    y = numeric(length(gms))
    for (i in seq_along(gms))
    {
      gm = gms[i]
      y[i] = integrate(integrand_sf,-Inf,Inf,gm,gf,baseline)$value
    }
    return(y)
  }
  integrand_gf = function(gfs,baseline)
  {
    y = numeric(length(gfs))
    for (i in seq_along(gfs))
    {
      gf = gfs[i]
      y[i] = integrate(integrand_gm,-Inf,Inf,gf,baseline)$value
    }
    return(y)
  }
  
  start_time = Sys.time()
  
  risk_selection = integrate(integrand_gf,-Inf,Inf,F)$value
  risk_baseline = integrate(integrand_gf,-Inf,Inf,T)$value
  
  end_time = Sys.time()
  cat(sprintf("Total time: %g\n",end_time-start_time))
  
  if (relative) {
    reduction = (risk_baseline-risk_selection)/risk_baseline
  } else {
    reduction = risk_baseline-risk_selection
  }
  return(reduction)
}