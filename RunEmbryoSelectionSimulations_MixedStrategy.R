
r2 = 0.1
Ks = c(0.01,0.05,0.2)

qs = seq(0,1,by=0.02)
n = 5

risk_red_highrisk_sim_all = array(0,c(length(Ks),length(qs),2))

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  cat(sprintf('K=%g\n',K))
  risk_red_highrisk_sim_all[ki,,1] = simulate_exclude_high(r2,K,n,qs,1000000,F,0,0,T,F)
  risk_red_highrisk_sim_all[ki,,2] = simulate_exclude_high(r2,K,n,qs,1000000,F,0,0,T,T)
}

saveRDS(risk_red_highrisk_sim_all, file = "./Data/high_risk_exclude_sim_mixed.rds")