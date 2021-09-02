
library(MASS)

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
      term2 = dnorm(u)*pnorm(beta,lower.tail=F)^(n-1) * internal_int
      y[i] = term1 + term2
    }
    return(y)
  }
  risk = integrate(integrand_u,-Inf,Inf)$value
  reduction = (K-risk)/K
  return(reduction)
}

# Tried to vectorize the integrands, but time saving was very minor
risk_reduction_exclude_vec = function(r2,K,q,n)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  integrand_t = function(t,u)
  {
    y = dnorm(t)*pnorm((zk-r/sqrt(2)*(u+t))/sqrt(1-r^2),lower.tail=F)
    return(y)
  }
  integrand_u = function(u)
  {
    beta = zq*sqrt(2)-u
    internal_int = integrate(integrand_t,-Inf,beta,u)$value
    denom = pnorm(beta)
    if (denom==0) {denom=1e-300} # Avoid dividing by zero
    numer = dnorm(u)*(1-pnorm(beta,lower.tail=F)^n) * internal_int
    term1 = numer/denom
    
    internal_int = integrate(integrand_t,beta,Inf,u)$value
    term2 = dnorm(u)*pnorm(beta,lower.tail=F)^(n-1) * internal_int
    y = term1 + term2
    return(y)
  }
  integrand_u_vec = Vectorize(integrand_u,"u")
  risk = integrate(integrand_u_vec,-Inf,Inf)$value
  reduction = (K-risk)/K
  return(reduction)
}

risk_reduction_exclude_conditional = function(r2,K,q,n,qf,qm,relative=T)
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
  term2 = pnorm(gamma,lower.tail=F)^(n-1) * internal_int

  risk = term1 + term2
  
  if (relative) {
    reduction = (baseline-risk)/baseline
  } else {
    reduction = baseline-risk
  }
  return(reduction)
}

risk_reduction_exclude_family_history = function(r2,h2,K,q,n,df,dm)
{
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  
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
  
  integrand_t = function(t,gm,gf)
  {
    arg = (zk-t*r/sqrt(2)-(gm+gf)/2)/sqrt(1-h^2/2-r^2/2)
    y = dnorm(t)*pnorm(arg,lower.tail=F)
    return(y)
  }
  
  integrand_c = function(cs,gm,gf)
  {
    y = numeric(length(cs))
    for (i in seq_along(cs))
    {
      c = cs[i]
      
      gamma = zq*sqrt(2) - c/(r/sqrt(2))
      internal_int = integrate(integrand_t,-Inf,gamma,gm,gf,rel.tol = 1e-6)$value
      denom = pnorm(gamma)
      if (denom==0) {denom=1e-300} # Avoid dividing by zero
      numer = (1-pnorm(gamma,lower.tail=F)^n) * internal_int
      term1 = numer/denom
      
      internal_int = integrate(integrand_t,gamma,Inf,gm,gf)$value
      term2 = pnorm(gamma,lower.tail=F)^(n-1) * internal_int
      
      y[i] = term1 + term2
      
      fc = dnorm(c,mean=r^2/h^2 * (gm+gf)/2, sd=r/h * sqrt((h^2-r^2)/2))
      y[i] = y[i] * fc
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
        y[i] = integrate(integrand_c,-Inf,Inf,gm,gf)$value
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
  
  risk_baseline = integrate(integrand_gf,-Inf,Inf,T)$value
  risk_selection = integrate(integrand_gf,-Inf,Inf,F)$value
  
  end_time = Sys.time()
  cat(sprintf("Total time: %g\n",end_time-start_time))
  
  relative_reduction = (risk_baseline-risk_selection)/risk_baseline
  abs_reduction = risk_baseline-risk_selection
  return(c(risk_baseline,risk_selection,relative_reduction,abs_reduction))
}
