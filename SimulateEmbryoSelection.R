
r2s = c(0.05,0.1,0.3)
Ks = c(0.01,0.05,0.2)

qs = seq(0,0.4,by=0.01)
n = 10
ns = seq(1,20)

risk_red_highrisk_sim_all = array(numeric(length(Ks)*length(r2s)*length(qs)),c(length(Ks),length(r2s),length(qs)))
risk_red_lowestrisk_sim_all = array(numeric(length(Ks)*length(r2s)*length(ns)),c(length(Ks),length(r2s),length(ns)))

for (ki in seq_along(Ks))
{
  K = Ks[ki]
  cat(sprintf('K=%g',K))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    cat(sprintf('\nr2=%g\n',r2))
    risk_red_highrisk_sim_all[ki,ri,] = simulate_exclude_high(r2,K,n,qs,1000000)
    # risk_red_lowestrisk_sim_all[ki,ri,] = simulate_lowest_risk(r2,K,ns,1000000)
  }
}

saveRDS(risk_red_highrisk_sim_all, file = "./Data/high_risk_exclude_sim_n10.rds")
# saveRDS(risk_red_lowestrisk_sim_all, file = "./Data/select_lowest_risk_sim.rds")