
r2s = c(0.05,0.1,0.3)
Ks = c(0.01,0.05,0.2)

n=5
num_cs_bins = 10000

out_name = "./Data/select_lowest_risk_personal_rrr_dist.rds"

risk_red_lowestrisk_th_all = array(0,c(length(Ks),length(r2s),num_cs_bins))

for (Ki in seq_along(Ks))
{
  K = Ks[Ki]
  cat(sprintf('K=%g\n',K))
  
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    cat(sprintf('r2=%g\n',r2))
    
    qs = seq(0,100,length.out = num_cs_bins+1)
    dq = qs[2]-qs[1]
    qs = qs[1:(length(qs)-1)]+dq/2
    
    for (qi in seq_along(qs))
    {
      q = qs[qi]
      # cat(sprintf('q=%g\n',q))
      # We compute the relative risk reduction by assuming both parents have the same score percentile q, because either way only the parental average matters
      rrr = risk_reduction_lowest_conditional(r2,K,n,q,q,parental_avg_given=T)
      risk_red_lowestrisk_th_all[Ki,ri,qi] = rrr
    }
  }
}

saveRDS(risk_red_lowestrisk_th_all, file = out_name)