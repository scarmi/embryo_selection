
source("EmbryoSelectionAnalysis_Exclude.R")
source("EmbryoSelectionAnalysis_Lowest.R")

r2s = c(0.05,0.1,0.3)
K = 0.05
h2 = 0.4

dfs = c(F,T,T)
dms = c(F,F,T)

qs = seq(0,0.4,by=0.01)
nh = 5
ns = seq(1,20)

out_name_sim_hre = "./Data/high_risk_exclude_sim_cond_parents.rds"
out_name_sim_lrp = "./Data/select_lowest_risk_sim_cond_parents.rds"
out_name_th_hre = "./Data/high_risk_exclude_th_cond_parents.rds"
out_name_th_lrp = "./Data/select_lowest_risk_th_cond_parents.rds"

risk_red_highrisk_sim_all = array(0,c(length(r2s),length(qs),length(dfs),2))
risk_red_lowestrisk_sim_all = array(0,c(length(r2s),length(ns),length(dfs),2))
risk_red_highrisk_th_all = array(0,c(length(r2s),length(qs),length(dfs),2))
risk_red_lowestrisk_th_all = array(0,c(length(r2s),length(ns),length(dfs),2))

for (ri in seq_along(r2s))
{
  r2 = r2s[ri]
  cat(sprintf('r2=%g\n',r2))
  print("Simulate exclude")
  res = simulate_exclude_high_parents_disease(r2,h2,K,nh,qs,1000000)
  
  risk_red_highrisk_sim_all[ri,,,1] = res[,,1]
  risk_red_highrisk_sim_all[ri,,,2] = res[,,2]
  
  print("Simulate lowest")
  res = simulate_lowest_risk_parents_disease(r2,h2,K,ns,1000000)
  risk_red_lowestrisk_sim_all[ri,,,1] = res[,,1]
  risk_red_lowestrisk_sim_all[ri,,,2] = res[,,2]
  
  for (di in seq_along(dfs))
  {
    df = dfs[di]
    dm = dms[di]
    cat(sprintf('df=%d, dm=%d\n',df,dm))
    for (qi in seq_along(qs))
    {
      q = qs[qi]
      cat(sprintf("theory exclude, q=%g\n",q))
      res = risk_reduction_exclude_family_history(r2,h2,K,q,nh,df,dm)
      risk_red_highrisk_th_all[ri,qi,di,1] = res[3]
      risk_red_highrisk_th_all[ri,qi,di,2] = res[4]
    }
    for (ni in seq_along(ns))
    {
      n = ns[ni]
      cat(sprintf("theory lowest, n=%g\n",n))
      res = risk_reduction_lowest_family_history(r2,h2,K,n,df,dm)
      risk_red_lowestrisk_th_all[ri,ni,di,1] = res[3]
      risk_red_lowestrisk_th_all[ri,ni,di,2] = res[4]
    }
  }
}

saveRDS(risk_red_highrisk_sim_all, file = out_name_sim_hre)
saveRDS(risk_red_lowestrisk_sim_all, file = out_name_sim_lrp)
saveRDS(risk_red_highrisk_th_all, file = out_name_th_hre)
saveRDS(risk_red_lowestrisk_th_all, file = out_name_th_lrp)