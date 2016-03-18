
calculate_summary_statistics <- function(data,pop1,pop2){
  # Heterozygosity
  H_pop1    <- apply(data[pop1,], 2, function(x) H(as.factor(x)) ) 
  m_H_pop1  <- mean(H_pop1,na.rm=T)
  v_H_pop1  <- var(H_pop1,na.rm=T)
  H_pop2    <- apply(data[pop2,], 2, function(x) H(as.factor(x)) ) 
  m_H_pop2  <- mean(H_pop2,na.rm=T)
  v_H_pop2  <- var(H_pop2,na.rm=T)
  H_total   <- apply(data, 2, function(x) H(as.factor(x)) ) 
  m_H_total <- mean(H_total,na.rm=T)
  v_H_total <- var(H_total,na.rm=T)
  # Number of alleles
  A_pop1    <- apply(data[pop1,], 2, function(x) length(levels(as.factor(x))) ) 
  m_A_pop1  <- mean(A_pop1,na.rm=T)
  v_A_pop1  <- var(A_pop1,na.rm=T)
  A_pop2    <- apply(data[pop2,], 2, function(x) length(levels(as.factor(x))) ) 
  m_A_pop2  <- mean(A_pop2,na.rm=T)
  v_A_pop2  <- var(A_pop2,na.rm=T)
  A_total   <- apply(data, 2, function(x) length(levels(as.factor(x))) ) 
  m_A_total <- mean(A_total,na.rm=T)
  v_A_total <- var(A_total,na.rm=T)
  # Number of private alleles
  PA_pop1    <- apply(data[,], 2, function(x) length(which(!is.na(setdiff(unique(x[pop1]),unique(x[pop2]))))) ) 
  m_PA_pop1  <- mean(PA_pop1,na.rm=T)
  v_PA_pop1  <- var(PA_pop1,na.rm=T)
  PA_pop2    <- apply(data[,], 2, function(x) length(which(!is.na(setdiff(unique(x[pop2]),unique(x[pop1]))))) ) 
  m_PA_pop2  <- mean(PA_pop2,na.rm=T)
  v_PA_pop2  <- var(PA_pop2,na.rm=T)
  # Number of shared alleles
  SA         <- apply(data[,], 2, function(x) length(which(!is.na(intersect(unique(x[pop1]),unique(x[pop2]))))) ) 
  m_SA       <- mean(SA,na.rm=T)
  v_SA       <- var(SA,na.rm=T)
  # Variance in allele size
  V_pop1    <- apply(data[pop1,], 2, var, na.rm = TRUE) 
  m_V_pop1  <- mean(V_pop1,na.rm=T)
  v_V_pop1  <- var(V_pop1,na.rm=T)
  V_pop2    <- apply(data[pop2,], 2, var, na.rm = TRUE) 
  m_V_pop2  <- mean(V_pop2,na.rm=T)
  v_V_pop2  <- var(V_pop2,na.rm=T)
  V_total   <- apply(data, 2, var, na.rm = TRUE) 
  m_V_total <- mean(V_total,na.rm=T)
  v_V_total <- var(V_total,na.rm=T)
  # Range in allele size
  R_pop1    <- apply(data[pop1,], 2, function(x) max(x,na.rm=T)-min(x,na.rm=T) ) 
  m_R_pop1  <- mean(R_pop1,na.rm=T)
  v_R_pop1  <- var(R_pop1,na.rm=T)
  R_pop2    <- apply(data[pop2,], 2, function(x) max(x,na.rm=T)-min(x,na.rm=T)) 
  m_R_pop2  <- mean(R_pop2,na.rm=T)
  v_R_pop2  <- var(R_pop2,na.rm=T)
  R_total   <- apply(data, 2, function(x) max(x,na.rm=T)-min(x,na.rm=T)) 
  m_R_total <- mean(R_total,na.rm=T)
  v_R_total <- var(R_total,na.rm=T)
  # Garza and Williamson statistic
  GW_pop1    <- A_pop1/(R_pop1+1)
  m_GW_pop1  <- mean(GW_pop1,na.rm=T)
  v_GW_pop1  <- var(GW_pop1,na.rm=T)
  GW_pop2    <- A_pop2/(R_pop2+1)
  m_GW_pop2  <- mean(GW_pop2,na.rm=T)
  v_GW_pop2  <- var(GW_pop2,na.rm=T)
  GW_total   <- A_total/(R_total+1)
  m_GW_total <- mean(GW_total,na.rm=T)
  v_GW_total <- var(GW_total,na.rm=T)
  # Ratio between heterozygosity and allele size variance
  HV_pop1    <- H_pop1/V_pop1
  m_HV_pop1  <- mean(HV_pop1,na.rm=T)
  v_HV_pop1  <- var(HV_pop1,na.rm=T)
  HV_pop2    <- H_pop2/V_pop2
  m_HV_pop2  <- mean(HV_pop2,na.rm=T)
  v_HV_pop2  <- var(HV_pop2,na.rm=T)
  HV_total   <- H_total/V_total
  m_HV_total <- mean(HV_total,na.rm=T)
  v_HV_total <- var(HV_total,na.rm=T)
  # GST
  GST_by_locus <- (H_total - apply(rbind(H_pop1,H_pop2),2,mean))/(H_total) 
  GST          <- sum(H_total - apply(rbind(H_pop1,H_pop2),2,mean))/sum(H_total)
  m_GST        <- mean(GST_by_locus,na.rm=T)
  v_GST        <- var(GST_by_locus,na.rm=T)
  # delta mu squared statistic
  du2   <- apply(data, 2, function(x) (mean(x[pop1],na.rm=T)-mean(x[pop2],na.rm=T))^2 ) 
  m_du2 <- mean(du2,na.rm=T)
  v_du2 <- var(du2,na.rm=T)
  
  return(cbind(m_H_pop1,  v_H_pop1,  m_H_pop2,  v_H_pop2,  m_H_total,  v_H_total,
               m_A_pop1,  v_A_pop1,  m_A_pop2,  v_A_pop2,  m_A_total,  v_A_total,
               m_PA_pop1, v_PA_pop1, m_PA_pop2, v_PA_pop2, m_SA,       v_SA,
               m_V_pop1,  v_V_pop1,  m_V_pop2,  v_V_pop2,  m_V_total,  v_V_total,
               m_R_pop1,  v_R_pop1,  m_R_pop2,  v_R_pop2,  m_R_total,  v_R_total,
               m_GW_pop1, v_GW_pop1, m_GW_pop2, v_GW_pop2, m_GW_total, v_GW_total,
               m_HV_pop1, v_HV_pop1, m_HV_pop2, v_HV_pop2, m_HV_total, v_HV_total,
               GST,m_GST, v_GST,m_du2,v_du2))
}