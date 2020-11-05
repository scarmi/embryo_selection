
relative = T
if (relative) 
{
  ymax=99
  text_pos=92
  ylabel = "Risk reduction (%)"
  leg_pos_hre = 'topright'
  leg_pos_lrp = 'bottomright'
  input_name_hre = "./Data/high_risk_exclude_sim_cond.rds"
  output_name_hre = "../Figures/ExcludeHigh_Conditional.png"
  input_name_lrp = './Data/select_lowest_risk_sim_cond.rds'
  output_name_lrp = "../Figures/LowestRisk_Conditional.png"
  
} else {
  ymax=28
  text_pos=26
  ylabel = "Risk reduction (% points)"
  leg_pos_hre = 'right'
  leg_pos_lrp = 'right'
  input_name_hre = "./Data/high_risk_exclude_sim_cond_abs.rds"
  output_name_hre = "../Figures/ExcludeHigh_Conditional_Absolute.png"
  input_name_lrp = './Data/select_lowest_risk_sim_cond_abs.rds'
  output_name_lrp = "../Figures/LowestRisk_Conditional_Absolute.png"
}

r2s = c(0.05,0.1,0.3)
K = 0.05
qfs = c(50,98,98,98)
qms = c(50,25,50,98)

colors = c('Maroon','DarkCyan','DarkSlateBlue')

baseline_risk = function(r2,K,qf,qm)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf/100)
  zqm = qnorm(qm/100)
  c = (zqf+zqm)/2 * r
  baseline= pnorm((zk-c)/sqrt(1-r^2/2),lower.tail=F)
  return(baseline)
}

######### Exclude high risk ############

qs = seq(0,0.4,by=0.01)
n = 5
labels = c('A','B','C','D')

risk_red_highrisk_sim_all = readRDS(file = input_name_hre)
png(output_name_hre,3200,1000,res=300)

par(mfrow=c(1,length(qfs)),mar=c(4,5.5,2,2),yaxs='i',xaxs='i')

for (qi in seq_along(qfs))
{
  qf = qfs[qi]
  qm = qms[qi]
  plot(0,type='n',xlab="Percentile PRS to exclude",ylab=ylabel,cex.lab=1.5,cex.axis=1.25,cex.main=1.15,main=sprintf('HR exclusion, par. 1=%g%%, par. 2=%g%%',qf,qm),font.main=1,ylim=c(0,ymax),xlim=c(0,max(qs)*100))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    risk_red_sim = risk_red_highrisk_sim_all[qi,ri,]
    risk_red_th = numeric(length(qs))
    for (i in seq_along(qs))
    {
      q = qs[i]
      risk_red_th[i] = risk_reduction_exclude_conditional(r2,K,q,n,qf,qm,relative)
    }
    points(qs*100,risk_red_sim*100,cex=1.3,lwd=1.5,col=colors[ri],pch=ri-1)
    lines(qs*100,risk_red_th*100,lwd=1.8,col=colors[ri],lty=1)
    
    index = 3
    segments(qs[index]*100,-10,qs[index]*100,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
    segments(qs[index]*100,risk_red_sim[index]*100,-10,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
  }
  legs = numeric(length(r2s))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    baseline = round(baseline_risk(r2,K,qf,qm),2)*100
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)*', bl='*.(baseline)*'%'))
    # legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }

  legend(leg_pos_hre,legs,col=colors,pch=0:2,bty='n',cex=1.2,text.col='DarkBlue',lwd=2)
  text(4,text_pos,labels[qi],cex=2)
}
dev.off()

######### Select lowest risk ############

ns = seq(1,20)
labels = c('E','F','G','H')

risk_red_lowestrisk_sim_all = readRDS(file = input_name_lrp)
png(output_name_lrp,3200,1000,res=300)

par(mfrow=c(1,length(qfs)),mar=c(4,5.5,2,2),yaxs='i',xaxs='i')

for (qi in seq_along(qfs))
{
  qf = qfs[qi]
  qm = qms[qi]
  plot(0,type="n",xlab="Number of embryos",ylab=ylabel,cex.lab=1.5,cex.axis=1.25,main=sprintf('Lowest risk, par. 1=%g%%, par. 2=%d%%',qf,qm),font.main=1,cex.main=1.15,ylim=c(0,ymax),xlim=c(1,20))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    risk_red_sim = risk_red_lowestrisk_sim_all[qi,ri,]
    risk_red_th = numeric(length(ns))
    for (i in seq_along(ns))
    {
      n = ns[i]
      risk_red_th[i] = risk_reduction_lowest_conditional(r2,K,n,qf,qm,relative)
    }
    points(ns,risk_red_sim*100,cex=1.3,lwd=1.5,col=colors[ri],pch=ri-1)
    lines(ns,risk_red_th*100,lwd=1.8,lty=1,col=colors[ri])
    
    index = 5
    segments(ns[index],-10,ns[index],risk_red_sim[index]*100,col='black',lty=1,lwd=2)
    segments(ns[index],risk_red_sim[index]*100,-10,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
  }
  legs = numeric(length(r2s))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    baseline = round(baseline_risk(r2,K,qf,qm),2)*100
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)*', bl='*.(baseline)*'%'))
    # legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }
  legend(leg_pos_lrp,legs,col=colors,pch=0:2,bty='n',cex=1.3,text.col='DarkBlue',lwd=2)
  text(3,text_pos,labels[qi],cex=2)
}
dev.off()
