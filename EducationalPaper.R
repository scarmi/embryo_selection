
for (r2 in c(0.05,0.1))
{
  for (K in c(0.005,0.01,0.1))
  {
    for (n in c(2,3,5))
    {
      rrr = risk_reduction_lowest(r2,K,n)
      nns = 1/(K*rrr)
      cat(sprintf("r2=%g, K=%g, n=%d. RRR=%.2g, NNS=%d\n",r2,K,n,rrr,round(nns)))
    }
  }
}