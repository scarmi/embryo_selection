
library(MASS)

simulate_lowest_risk_two_traits = function(r2A,r2B,rho,KA,KB,ns,nfam=10000)
{
  results = array(0,c(length(ns),4))
  
  disease_countA = numeric(length(ns))
  disease_count_randomA = numeric(length(ns))
  disease_countB = numeric(length(ns))
  disease_count_randomB = numeric(length(ns))
  
  baselineA = numeric(length(ns))
  baselineB = numeric(length(ns))
  
  t_sickA = qnorm(1-KA)
  t_sickB = qnorm(1-KA)
  
  rA = sqrt(r2A)
  rB = sqrt(r2B)
  
  for (i in seq_along(ns))
  {
    cat('\r',i)
    
    n = ns[i]
    disease_countA[i] = 0
    disease_countB[i] = 0
    disease_count_randomA[i] = 0
    disease_count_randomB[i] = 0
    
    mu = c(0,0)
    sigma = 1/2 * matrix(c(r2A,rho*rA*rB,rho*rA*rB,r2B),nrow=2)
    
    cs = mvrnorm(nfam,mu,sigma)
    csA = cs[,1]
    csB = cs[,2]
    
    xs = mvrnorm(nfam*n,mu,sigma)
    xsA = matrix(xs[,1],nrow=nfam)
    xsB = matrix(xs[,2],nrow=nfam)
    
    envsA = rnorm(nfam,0,sqrt(1-r2A))
    envsB = rnorm(nfam,0,sqrt(1-r2B))
    
    for (j in seq(1,nfam))
    {
      cA = csA[j]
      cB = csB[j]
      
      xA = xsA[j,]
      xB = xsB[j,]
      
      scoresA = xA+cA
      scoresB = xB+cB
      
      envA = envsA[j]
      envB = envsB[j]
      
      selected_ind = which.min(scoresA)
      
      liabA = scoresA[selected_ind] + envA
      liabB = scoresB[selected_ind] + envB
      
      if (liabA>t_sickA)
        disease_countA[i] = disease_countA[i]+1
      if (liabB>t_sickB)
        disease_countB[i] = disease_countB[i]+1
      
      # Liability when not selecting
      liab_randomA = scoresA[1]+envA
      if (liab_randomA>t_sickA)
        disease_count_randomA[i] = disease_count_randomA[i]+1
      liab_randomB = scoresB[1]+envB
      if (liab_randomB>t_sickB)
        disease_count_randomB[i] = disease_count_randomB[i]+1
    }

    rrrA = (disease_count_randomA[i] - disease_countA[i])/disease_count_randomA[i]
    arrA = (disease_count_randomA[i] - disease_countA[i])/nfam
    rrrB = (disease_count_randomB[i] - disease_countB[i])/disease_count_randomB[i]
    arrB = (disease_count_randomB[i] - disease_countB[i])/nfam
    
    results[i,1] = rrrA
    results[i,2] = arrA
    results[i,3] = rrrB
    results[i,4] = arrB
  }  
  return(results)
}

r2 = 0.1
Ks = c(0.01,0.05,0.2)
rhos = c(-0.3,seq(-0.2,-0.05,by=0.05))

ns = seq(1,20)

risk_red_lowest_sim_all = array(0,c(length(Ks),length(ns),length(rhos)+1))

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  cat(sprintf('\nK=%g\n',K))
  for (rhi in seq_along(rhos))
  {
    rho = rhos[rhi]
    cat(sprintf('\nrho=%g\n',rho))
    res = simulate_lowest_risk_two_traits(r2,r2,rho,K,K,ns,1e6)
    if (rhi==1)
    {
      risk_red_lowest_sim_all[ki,,1] = res[,1]
    }
    risk_red_lowest_sim_all[ki,,rhi+1] = res[,3]
  }
}

saveRDS(risk_red_lowest_sim_all, file = "./Data/lowest_risk_two_taits.rds")


