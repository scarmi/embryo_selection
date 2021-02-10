
library(MASS)

simulate_exclude_high = function(r2,K,n,qs,nfam=10000,parents_known,qf,qm,relative=T,mixed=F)
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
        # If at least one embryo has score under the cutoff
        if (any(scores<t_exclude))
        {
          s = min(which(scores<t_exclude)) # Take the first embryo with score under the cutoff
        } else {
          # If all embryos have scores above the cutoff
          if (mixed) # Mixed strategy: in this case take the embryo with the lowest score 
          {
            s = which.min(scores)
          } else { # Otherwise take the first embryo
            s = 1  
          }
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