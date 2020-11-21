
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
    #print(t)
    #cat(sprintf("gm=%g, gf=%g\n",gm,gf))
    arg = (zk-t*r/sqrt(2)-(gm+gf)/2)/sqrt(1-h^2/2-r^2/2)
    #print(arg)
    y = dnorm(t)*pnorm(arg,lower.tail=F)
    #print(y)
    return(y)
  }
  
  integrand_c = function(cs,gm,gf)
  {
    y = numeric(length(cs))
    for (i in seq_along(cs))
    {
      c = cs[i]
      
      gamma = zq*sqrt(2) - c/(r/sqrt(2))
      # cat(sprintf("gamma=%g, c=%g, gm=%g, gf=%g\n",gamma,c,gm,gf))
      internal_int = integrate(integrand_t,-Inf,gamma,gm,gf,rel.tol = 1e-5)$value
      #cat(sprintf("internal_int1=%g\n",internal_int))
      denom = pnorm(gamma)
      if (denom==0) {denom=1e-300} # Avoid dividing by zero
      numer = (1-pnorm(gamma,lower.tail=F)^n) * internal_int
      term1 = numer/denom
      
      internal_int = integrate(integrand_t,gamma,Inf,gm,gf)$value
      #cat(sprintf("internal_int2=%g\n",internal_int))
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

simulate_exclude_high = function(r2,K,n,qs,nfam=10000,parents_known,qf,qm,relative=T)
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
  if (relative) {
    risk_red_sim = (baseline-disease_count/nfam)/baseline
  } else {
    risk_red_sim = baseline-disease_count/nfam
  }
  return(risk_red_sim)
}

simulate_exclude_high_parents_disease = function(r2,h2,K,n,qs,nfam=10000)
{
  num_histories = 3
  scores = numeric(nfam)
  results = array(0,c(length(qs),num_histories,4))
  
  t_sick = qnorm(1-K)
  r = sqrt(r2)
  h = sqrt(h2)
  
  for (qi in seq_along(qs))
  {
    cat('\r',qi)
    
    q = qs[qi]
    t_exclude = qnorm(1-q)*sqrt(r2)
    
    disease_count = numeric(num_histories)
    disease_count_random = numeric(num_histories)
    family_count = numeric(num_histories)
    
    sm = rnorm(nfam,0,r)
    gtm = rnorm(nfam,0,sqrt(h^2-r^2))
    sf = rnorm(nfam,0,r)
    gtf = rnorm(nfam,0,sqrt(h^2-r^2))
    envm = rnorm(nfam,0,sqrt(1-h^2))
    envf = rnorm(nfam,0,sqrt(1-h^2))
    
    c = (sm+sf)/2
    cg = (gtm+gtf)/2
    
    xs = rnorm(nfam*n,0,r/sqrt(2))
    xs = matrix(xs,nrow=nfam)
    gtc = rnorm(nfam,0,sqrt((h^2-r^2)/2))
    envc = rnorm(nfam,0,sqrt(1-h^2))
    
    for (j in seq(1,nfam))
    {
      liabm = sm[j] + gtm[j] + envm[j]
      liabf = sf[j] + gtf[j] + envf[j]
      count_sick = 0
      if (liabm >= t_sick) {count_sick = count_sick + 1}
      if (liabf >= t_sick) {count_sick = count_sick + 1}
      history = count_sick + 1 # 1: both healthy, 2: one sick, 3: both sick
      family_count[history] = family_count[history]+1
      
      if (n==Inf)
      {
        while (T)
        {
          x = rnorm(1,0,sqrt(r2/2))
          score = x+c[j]
          if (score<t_exclude)
            break
        }
        liabc = score + cg[j] + gtc[j] + envc[j]
      } else {
        x = xs[j,]
        scores = x+c[j]
        if (any(scores<t_exclude))
        {
          s = min(which(scores<t_exclude))
        } else {
          s = 1
        }
        liabc = scores[s] + cg[j] + gtc[j] + envc[j]
      }
      if (liabc>t_sick)
      {
        disease_count[history] = disease_count[history]+1
      }
      # Liability when not selecting
      x = rnorm(1,0,sqrt(r2/2))
      liabc = c[j] + x[1] + cg[j] + gtc[j] + envc[j]
      if (liabc>t_sick)
      {
        disease_count_random[history] = disease_count_random[history]+1
      }
    }
    risk_selection = disease_count / family_count
    risk_random = disease_count_random / family_count
    rrr = (risk_random - risk_selection)/risk_random
    arr = (risk_random - risk_selection)
    results[qi,,1] = rrr
    results[qi,,2] = arr
    results[qi,,3] = risk_selection
    results[qi,,4] = risk_random
  }
  return(results)
}