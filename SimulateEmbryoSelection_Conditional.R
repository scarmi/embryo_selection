
r2s = c(0.05,0.1,0.3)
K = 0.05

qfs = c(50,98,98,98)
qms = c(50,25,50,98)

qs = seq(0,0.4,by=0.01)
n = 5
ns = seq(1,20)

relative=F

if (relative) {
  out_name_hre = "./Data/high_risk_exclude_sim_cond.rds"
  out_name_lrp = "./Data/select_lowest_risk_sim_cond.rds"
} else {
  out_name_hre = "./Data/high_risk_exclude_sim_cond_abs.rds"
  out_name_lrp = "./Data/select_lowest_risk_sim_cond_abs.rds"
}

risk_red_highrisk_sim_all = array(numeric(length(qfs)*length(r2s)*length(qs)),c(length(qfs),length(r2s),length(qs)))
risk_red_lowestrisk_sim_all = array(numeric(length(qfs)*length(r2s)*length(ns)),c(length(qfs),length(r2s),length(ns)))

for (qi in seq_along(qfs))
{
  qf = qfs[qi]
  qm = qms[qi]
  
  cat(sprintf('qf=%g, qm=%g',qf,qm))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    cat(sprintf('\nr2=%g\n',r2))
    risk_red_highrisk_sim_all[qi,ri,] = simulate_exclude_high(r2,K,n,qs,1000000,T,qf,qm,relative)
    risk_red_lowestrisk_sim_all[qi,ri,] = simulate_lowest_risk(r2,K,ns,1000000,T,qf,qm,relative)
  }
}

saveRDS(risk_red_highrisk_sim_all, file = out_name_hre)
saveRDS(risk_red_lowestrisk_sim_all, file = out_name_lrp)