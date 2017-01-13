

# functions for additional summary statistics
Kimmel_beta <- function(He,V,log_beta=F){
  stopifnot(length(V)==length(He))
  P0 <- 1-He
  theta_V <- V
  theta_P <- (1/P0^2-1)/2
  beta <- theta_V/theta_P
  if (log_beta){
    beta <- log(beta)
  }
  return(beta)
}
diffNa <- function(M,n,obsNa){
  j <- 1:n-1
  expNa <- sum(M/(M+j))
  diffNa <- obsNa-expNa
  return(abs(diffNa))
}
#data_sample_size
bottleneck <- function(Na,He,n){
  M <- optimize(f=diffNa, interval=c(10^-6,10^6),n=n,obsNa=Na)$minimum
  expHe <- M/(1+M)
  bottleneck <- He - expHe
  return(bottleneck)
}
Navascues_beta <- function(Na,V,log_beta=F){
  stopifnot(length(V)==length(Na))
  theta_V  <- 2*V
  theta_Na <- (Na^2-1)/2
  beta <- theta_V/theta_Na
  if (log_beta){
    beta <- log(beta)
  }
  return(beta)
}







calculate_summary_statistics <- function(data,pop1,pop2){
  
  # Change to other formats
  data_as_loci <- data.frame(apply(data, 2, function(x) paste(x[seq(1,680,2)],x[seq(2,680,2)],sep="/")))
  data_as_loci[data_as_loci=="NA/NA"]<-NA
  data_as_genind <- df2genind(data_as_loci,sep="/",pop=c(rep(1,length(pop1)/2),rep(2,length(pop2)/2)))
  clusters<-c(rep(1,length(pop1)/2),rep(2,length(pop2)/2))
  data_as_loci <- cbind("population"=clusters,data_as_loci)
  data_as_loci <- as.loci(data_as_loci, allele.sep = "/")

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
  # Kimmel beta
  Beta_pop1   <- Kimmel_beta(H_pop1,V_pop1,log_beta=T)
  m_Beta_pop1 <- mean(Beta_pop1,na.rm=T)
  v_Beta_pop1 <- var(Beta_pop1,na.rm=T)
  Beta_pop2   <- Kimmel_beta(H_pop2,V_pop2,log_beta=T)
  m_Beta_pop2 <- mean(Beta_pop2,na.rm=T)
  v_Beta_pop2 <- var(Beta_pop2,na.rm=T)
  Beta_total   <- Kimmel_beta(H_total,V_total,log_beta=T)
  m_Beta_total <- mean(Beta_total,na.rm=T)
  v_Beta_total <- var(Beta_total,na.rm=T)
  # Navascués beta
  BetaP_pop1   <- Navascues_beta(A_pop1,V_pop1,log_beta=T)
  m_BetaP_pop1 <- mean(BetaP_pop1,na.rm=T)
  v_BetaP_pop1 <- var(BetaP_pop1,na.rm=T)
  BetaP_pop2   <- Navascues_beta(A_pop2,V_pop2,log_beta=T)
  m_BetaP_pop2 <- mean(BetaP_pop2,na.rm=T)
  v_BetaP_pop2 <- var(BetaP_pop2,na.rm=T)
  BetaP_total   <- Navascues_beta(A_total,V_total,log_beta=T)
  m_BetaP_total <- mean(BetaP_total,na.rm=T)
  v_BetaP_total <- var(BetaP_total,na.rm=T)
  # Bottleneck
  bottleneck_pop1  <- array(NA,ncol(data))
  bottleneck_pop2  <- array(NA,ncol(data))
  bottleneck_total <- array(NA,ncol(data))
  for (locus in seq_len(ncol(data)) ){
    popsize1 <- sum(!is.na(data[pop1,locus]))
    popsize2 <- sum(!is.na(data[pop2,locus]))
    popsizeT <- popsize1+popsize2
    bottleneck_pop1[locus]  <- bottleneck(A_pop1[locus],H_pop1[locus],popsize1)
    bottleneck_pop2[locus]  <- bottleneck(A_pop2[locus],H_pop2[locus],popsize2)
    bottleneck_total[locus] <- bottleneck(A_total[locus],H_total[locus],popsizeT)
  }
  m_bottleneck_pop1  <- mean(bottleneck_pop1,na.rm=T)
  v_bottleneck_pop1  <- var(bottleneck_pop1,na.rm=T)
  m_bottleneck_pop2  <- mean(bottleneck_pop2,na.rm=T)
  v_bottleneck_pop2  <- var(bottleneck_pop2,na.rm=T)
  m_bottleneck_total <- mean(bottleneck_total,na.rm=T)
  v_bottleneck_total <- var(bottleneck_total,na.rm=T)
  

  # delta mu squared statistic
  du2   <- apply(data, 2, function(x) (mean(x[pop1],na.rm=T)-mean(x[pop2],na.rm=T))^2 ) 
  m_du2 <- mean(du2,na.rm=T)
  v_du2 <- var(du2,na.rm=T)

  # Differentiation statistics as calculated by mmod package
  mmod_stats <- diff_stats(data_as_genind,phi_st=T)
  
  g_Gst       <- mmod_stats$global["Gst_est"]
  m_Gst       <- mean(mmod_stats$per.locus[,"Gst"],na.rm=T)
  v_Gst       <- var(mmod_stats$per.locus[,"Gst"],na.rm=T)

  g_Gst_prime <- mmod_stats$global["Gprime_st"]
  m_Gst_prime <- mean(mmod_stats$per.locus[,"Gprime_st"],na.rm=T)
  v_Gst_prime <- var(mmod_stats$per.locus[,"Gprime_st"],na.rm=T)
  
  g_JostD     <- mmod_stats$global["D_het"]
  h_JostD     <- mmod_stats$global["D_mean"]
  m_JostD     <- mean(mmod_stats$per.locus[,"D"],na.rm=T)
  v_JostD     <- var(mmod_stats$per.locus[,"D"],na.rm=T)

  g_Phi_st    <- mmod_stats$global["Phi_st"]
  m_Phi_st    <- mean(mmod_stats$per.locus[,"Phi_st"],na.rm=T)
  v_Phi_st    <- var(mmod_stats$per.locus[,"Phi_st"],na.rm=T)
  
  # Differentiation statistics from pegas package (Weir and Cockerham)
  F_stats    <- Fst(data_as_loci)

  m_Fst <- mean(F_stats[,"Fst"],na.rm=T)
  v_Fst <- var(F_stats[,"Fst"],na.rm=T)

  # Genetic distances from adegenet package
  Nei_dist      <- dist.genpop(genind2genpop(data_as_genind,quiet=T),method=1)
  Edwards_dist  <- dist.genpop(genind2genpop(data_as_genind,quiet=T),method=2)
  Reynolds_dist <- dist.genpop(genind2genpop(data_as_genind,quiet=T),method=3)
  Rogers_dist   <- dist.genpop(genind2genpop(data_as_genind,quiet=T),method=4)
  Provesti_dist <- dist.genpop(genind2genpop(data_as_genind,quiet=T),method=5)
  
  
  
  
  
    
  return(cbind(m_H_pop1,  v_H_pop1,  m_H_pop2,  v_H_pop2,  m_H_total,  v_H_total,
               m_A_pop1,  v_A_pop1,  m_A_pop2,  v_A_pop2,  m_A_total,  v_A_total,
               m_PA_pop1, v_PA_pop1, m_PA_pop2, v_PA_pop2, m_SA,       v_SA,
               m_V_pop1,  v_V_pop1,  m_V_pop2,  v_V_pop2,  m_V_total,  v_V_total,
               m_R_pop1,  v_R_pop1,  m_R_pop2,  v_R_pop2,  m_R_total,  v_R_total,
               m_GW_pop1, v_GW_pop1, m_GW_pop2, v_GW_pop2, m_GW_total, v_GW_total,
               m_Beta_pop1,  v_Beta_pop1,  m_Beta_pop2,  v_Beta_pop2,  m_Beta_total,  v_Beta_total,
               m_BetaP_pop1, v_BetaP_pop1, m_BetaP_pop2, v_BetaP_pop2, m_BetaP_total, v_BetaP_total,
               m_bottleneck_pop1, v_bottleneck_pop1,
               m_bottleneck_pop2, v_bottleneck_pop2,
               m_bottleneck_total, v_bottleneck_total,
               g_Gst,            m_Gst,       v_Gst,
               g_Gst_prime,      m_Gst_prime, v_Gst_prime,
               g_JostD, h_JostD, m_JostD,     v_JostD,
               g_Phi_st,         m_Phi_st,    v_Phi_st,
               Nei_dist, Edwards_dist, Reynolds_dist, Rogers_dist, Provesti_dist,
               m_Fst, v_Fst,
               m_du2, v_du2))
}


do_sim <- function(sim,nsim,
                   theta_min,theta_max,
                   M_min,M_max,
                   T_min,T_max,
                   alpha_min,alpha_max,
                   PGSM_min,PGSM_max,
                   sample_size_total,sample_size_pop1,sample_size_pop2,
                   num_of_loci, missing_data){
  # write progress of simulation on screen 
  if (sim==1 | sim%%10==0 | sim==nsim){
    cat(paste(sim,"of",nsim))  
  }
  
  # SAMPLE FROM PRIOR
  #-------------------
  
  # model (migration OR pure divergence) is sampled from prior
  model <- sample(c(0,1),1)
  
  # take parameter values from priors
  #theta_pop1 <- 10^runif( 2, min=theta_min, max=theta_max )
  #theta1 <- max(theta_pop1)
  #thetaF <- min(theta_pop1)
  theta1 <- 10^runif( 1, min=theta_min, max=theta_max )
  theta2 <- 10^runif( 1, min=theta_min, max=theta_max )
  #thetaA <- 10^runif( 1, min=theta_min, max=theta_max )
  x2     <- theta2/theta1
  #xA     <- thetaA/theta1
  #xF     <- thetaF/theta1
  #alpha1  <- runif( 1, min=alpha_min, max=alpha_max)
  #alpha2 <- runif( 1, min=alpha_min, max=alpha_max)
  if (model==1){
    M <- 10^runif( 1, min=M_min, max=M_max )
    #mig12   <- runif( 1, min=M_min, max=M_max )
    #mig21   <- runif( 1, min=M_min, max=M_max )
  }else{
    M <- 0
    #mig21<-mig12<-0
  }
  TS <- 10^runif( 1, min=T_min, max=T_max )
  #T1 <- runif( 1, min=T_min, max=TS )
  #T2 <- runif( 1, min=T_min, max=T_max )
  PGSM <- runif( 1, min=PGSM_min, max=PGSM_max )
  
  # RUM MS
  #--------
  
  # generate text with ms command
  ms_out_file <- paste0("msout",sim,".txt")
  ms_run <- paste( "bin/ms", sample_size_total, num_of_loci,          # total sample size & number of loci
                   "-t", theta1,                                  # population size in pop1
                   "-I 2", sample_size_pop1, sample_size_pop2, M, # sample sizes per population & migration
                   "-n 2", x2,                                    # population size in pop2
                   #"-g 1", alpha1,
                   #"-en", T1, "1", xF,                           # change pop size in pop1
                   #"-en", T2, "2", xA,                            # change pop size in pop2
                   "-ej", TS, "1 2",                              # creation of pop1 from pop2
                   ">", ms_out_file)                              # output file
  
  # rum ms on system
  if(.Platform$OS.type == "unix") {
    ms_run <- paste( "./", ms_run, sep="" )
    system( ms_run )
  }
  
  # READ SIMULATED DATA AND TRANSFORM TO MICROSATELLITE
  #-----------------------------------------------------
  
  # read simulation output (generated with an infinite site model, ISM)
  msout <- read.ms.output(file=ms_out_file)
  system( paste("rm",ms_out_file) )
  
  # transform simulated data from ISM into microsatellite evolving under GSM
  # (generalised stepwise mutation model)
  mutations <- msout$segsites
  if (sum(mutations)>0){
    simulated_data <- matrix(NA, nrow=sample_size_total, ncol=num_of_loci)
    for(locus in 1:num_of_loci){
      mutation_size <- rgeom(mutations[locus],1-PGSM)+1
      mutation_size <- mutation_size * sample(c(-1,1),mutations[locus],replace=T)
      simulated_data[,locus] <- 100+colSums(t(msout$gametes[[locus]])*mutation_size)
    }
    # add missing data to simulated data
    for(i in 1:nrow(missing_data)){
      simulated_data[missing_data[i,1],missing_data[i,2]]<-NA
    }
    remove(mutations)
    
    # CALCULATE SUMMARY STATISTICS FOR SIMULATED DATA
    #-------------------------------------------------
    
    sim_sumstats <- calculate_summary_statistics(data=simulated_data,
                                                 1:sample_size_pop1,
                                                 sample_size_pop1+1:sample_size_pop2)
    
    # WRITE PARAMATERS AND SUMMARY STATISTICS TO REFERENCE TABLE
    #------------------------------------------------------------
    return( cbind( ## PARAMETERS
                  model,      # 0:isolation; 1:isolation with migration
                  theta1,     # theta1 (pop1)
#                  thetaF,     # theta1 (pop1)
                  theta2,     # theta2 (pop2)
                  thetaA,     # thetaA (ancestral)
                  alpha1,     # growth rate in pop1
#                  T1,         # change pop size in pop1
                  T2,         # change pop size in pop2
                  M,          # migration
                  TS,         # time split of pop1 from pop2
                  PGSM,
                  #### SUMMARY STATISTICS                      
                  sim_sumstats))
  }
}
  
  
  