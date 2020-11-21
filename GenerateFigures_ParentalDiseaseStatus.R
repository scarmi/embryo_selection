
relative = T
if (relative) 
{
  ymax=99
  text_pos=92
  ylabel = "Risk reduction (%)"
  leg_pos_hre = 'topright'
  leg_pos_lrp = 'bottomright'
  output_name_hre = "../Figures/ExcludeHigh_ConditionalParentalStatus.png"
  output_name_lrp = "../Figures/LowestRisk_ConditionalParentalStatus.png"
  
} else {
  ymax=28
  text_pos=26
  ylabel = "Risk reduction (% points)"
  leg_pos_hre = 'right'
  leg_pos_lrp = 'right'
  output_name_hre = "../Figures/ExcludeHigh_ConditionalParentalStatus_Absolute.png"
  output_name_lrp = "../Figures/LowestRisk_ConditionalParentalStatus_Absolute.png"
}

input_name_sim_hre_prefix = "./Data/Parental/high_risk_exclude_sim_cond_parents"
input_name_th_hre_prefix = "./Data/Parental/high_risk_exclude_th_cond_parents"
input_name_sim_lrp = './Data/select_lowest_risk_sim_cond_parents.rds'
input_name_th_lrp = './Data/select_lowest_risk_th_cond_parents.rds'

r2s = c(0.05,0.1,0.3)
K = 0.05
h2 = 0.4
dms = c(F,T,T)
dfs = c(F,F,T)

desc = c('both Healthy','one sick','both sick')
colors = c('Maroon','DarkCyan','DarkSlateBlue')

######### Exclude high risk ############
  
qs = seq(0,0.4,by=0.02)
n = 5
labels = c('A','B','C')

png(output_name_hre,3200,1000,res=300)

par(mfrow=c(1,length(dfs)),mar=c(4,5.5,2,2),yaxs='i',xaxs='i')

for (di in seq_along(dfs))
{
  df = dfs[di]
  dm = dms[di]
  plot(0,type='n',xlab="Percentile PRS to exclude",ylab=ylabel,cex.lab=1.5,cex.axis=1.25,cex.main=1.3,main=sprintf('HR exclusion, %s',desc[di]),font.main=1,ylim=c(0,ymax),xlim=c(0,max(qs)*100))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    
    risk_red_sim = numeric(length(qs))
    risk_red_th = numeric(length(qs))
    for (qi in seq_along(qs))
    {
      q = qs[qi]
      
      risk_red_highrisk_sim = readRDS(file = sprintf("%s_r_%d_q_%d.rds",input_name_sim_hre_prefix,ri,qi))
      risk_red_highrisk_th = readRDS(file = sprintf("%s_r_%d_q_%d.rds",input_name_th_hre_prefix,ri,qi))
      if (relative) {
        risk_red_sim[qi] = risk_red_highrisk_sim[di,1]
        risk_red_th[qi] = risk_red_highrisk_th[di,1]
      } else {
        risk_red_sim[qi] = risk_red_highrisk_sim[di,2]
        risk_red_th[qi] = risk_red_highrisk_th[di,2]
      }
    }
    points(qs*100,risk_red_sim*100,cex=1.3,lwd=1.5,col=colors[ri],pch=ri-1)
    lines(qs*100,risk_red_th*100,lwd=1.8,col=colors[ri],lty=1)
    
    index = 2
    segments(qs[index]*100,-10,qs[index]*100,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
    segments(qs[index]*100,risk_red_sim[index]*100,-10,risk_red_sim[index]*100,col='black',lty=1,lwd=2)
  }
  legs = numeric(length(r2s))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }
  
  legend(leg_pos_hre,legs,col=colors,pch=0:2,bty='n',cex=1.2,text.col='DarkBlue',lwd=2)
  
  # Here we're using the lowest risk strategy, but just for computing the baseline
  baseline = risk_reduction_lowest_family_history(0,h2,K,1,df,dm)[1]*100
  text(11,text_pos,sprintf('%s, baseline=%.1f%%',labels[di],baseline),cex=1.3)
}
dev.off()

######### Select lowest risk ############

ns = seq(1,20)
labels = c('D','E','F')

risk_red_lowestrisk_sim_all = readRDS(file = input_name_sim_lrp)
risk_red_lowestrisk_th_all = readRDS(file = input_name_th_lrp)
png(output_name_lrp,3200,1000,res=300)

par(mfrow=c(1,length(dfs)),mar=c(4,5.5,2,2),yaxs='i',xaxs='i')

for (di in seq_along(dfs))
{
  df = dfs[di]
  dm = dms[di]
  plot(0,type="n",xlab="Number of embryos",ylab=ylabel,cex.lab=1.5,cex.axis=1.25,main=sprintf('Lowest risk, %s',desc[di]),font.main=1,cex.main=1.3,ylim=c(0,ymax),xlim=c(1,20))
  for (ri in seq_along(r2s))
  {
    r2 = r2s[ri]
    if (relative) {index=1} else {index=2}
    risk_red_sim = risk_red_lowestrisk_sim_all[ri,,di,index]
    risk_red_th = risk_red_lowestrisk_th_all[ri,,di,index]
  
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
    legs[ri] = as.expression(bquote('r'^'2'*'='*.(r2)))
  }
  legend(leg_pos_lrp,legs,col=colors,pch=0:2,bty='n',cex=1.3,text.col='DarkBlue',lwd=2)
  
  baseline = risk_reduction_lowest_family_history(0,h2,K,1,df,dm)[1]*100
  text(6,text_pos,sprintf('%s, baseline=%.1f%%',labels[di],baseline),cex=1.3)
}
dev.off()
