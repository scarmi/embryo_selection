
library(MASS)

simulate_lowest_risk = function(r2,K,ns,nfam=10000,parents_known=F,qf,qm,relative=T)
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
  if (relative) {
    risk_red_sim = (baseline-disease_count/nfam)/baseline
  } else {
    risk_red_sim = baseline-disease_count/nfam
  }
  return(risk_red_sim)
}

simulate_lowest_risk_parents_disease = function(r2,h2,K,ns,nfam=10000)
{
  num_histories = 3
  scores = numeric(nfam)
  results = array(0,c(length(ns),num_histories,4))
  
  t_sick = qnorm(1-K)
  r = sqrt(r2)
  h = sqrt(h2)
  
  for (ni in seq_along(ns))
  {
    cat('\r',ni)
    
    n = ns[ni]
    
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
      
      x = xs[j,]
      liabc = c[j] + min(x) + cg[j] + gtc[j] + envc[j]
      if (liabc>t_sick)
        disease_count[history] = disease_count[history]+1
      # Liability when not selecting
      liabc = c[j] + x[1] + cg[j] + gtc[j] + envc[j]
      if (liabc>t_sick)
        disease_count_random[history] = disease_count_random[history]+1
    }
    
    risk_selection = disease_count / family_count
    risk_random = disease_count_random / family_count
    rrr = (risk_random - risk_selection)/risk_random
    arr = (risk_random - risk_selection)
    results[ni,,1] = rrr
    results[ni,,2] = arr
    results[ni,,3] = risk_selection
    results[ni,,4] = risk_random
  }
  return(results)
}