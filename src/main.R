#!/usr/bin/env Rscript

# Script to perform an ABC analysis
# Miguel Navascués 2017
# Miguel.Navascues@supagro.inra.fr

### Remove last features
rm(list=ls())
ls()


#########################################
# ABC GENERAL SETTINGS
#########################################

nsim             <- 200000               # number of simulations
batch_size       <- 100
nsim <- round(nsim/batch_size)*batch_size
input_filename   <- "data/agrarius.str" # data file
reftable_file    <- "results/reference_table"   # reference table file name
seed             <- 1234 ; 
set.seed(seed,"Mersenne-Twister")
check_data       <- T
simulations_only <- F
parallel_sims    <- T
num_of_threads   <- 10

# Model graphical description
#
#     pop1            pop2
#
#    theta1          theta2
#     |   |    M    |     | 
#     |   |  <--->  |     | 
#     |   |         |     | 
#     |   |         |     |
#  TS |___|_________|_____|
#            |  |
#            |  |
#            |  |
#            |  |
# 
#
# NB: a population expansion is assumed in pop 1 (with growth rate alpha1),
# but no assumption on the type of demographic change
# is made for pop 2 (see file fun.R to check ms command)
# order of events (TS and T2) is not fixed a priori

#########################################
# DEFINE PRIORS (upper and lower limit)
#########################################

theta_min <- -1
theta_max <-  2

M_min <- -5
M_max <- 3

T_min <- -5
T_max <- 4

alpha_min <- 0
alpha_max <- 5

PGSM_min <- 0     # Generalised Stepwise Mutation model parameter 
PGSM_max <- 1     # (disstribution probability for size of mutations in repeat units)

#########################################
# LOAD REQUIRED FUNCTIONS/PACKAGES
#########################################

# load required libraries
library(plyr)   # for function count() that helps to explore allele size distribution
library(pegas)  # for calculating heterozygosity: function H()
library(mmod)   # for Jost's D, etc
library(abc)    # for approximate Bayesian computation functions
require(hexbin) # for plots of PCA
require(grid)   # for plots of PCA
require(abcrf)           # abc random forest: model choice
require(quantregForest)  # abc random forest: parameter estimation
if(parallel_sims){
  library(foreach)    #
  library(parallel)   # -> for simulating in parallel
  library(doParallel) #
}

# Uses RR Hudson coalescent simulator (ms)
# load RR Hudson's R script distributed with ms to read ms output in R
# http://home.uchicago.edu/~rhudson1/source/mksamples.html
source("bin/readms.output.R")

# load other functions (from file distributed together with present file)
source("src/fun.R")

#########################################
# LOAD DATA
#########################################

# read real data file
# file format: two-row STRUCTURE (*.str)
# first column: individual ID (two rows per diploid individual)
# second column: population ID (NB: only two populations should be in the input!)
# third...final column: alleles for locus 1...locus n
# alleles are coded as number of repeat motifs (NB: do not use size in bp!)
microsat_data <- read.table(input_filename,header=T)

# deal with missing data (change -9 to NA; record NA positions in data table)
head(microsat_data)
for(locus in 3:ncol(microsat_data)){
  microsat_data[which(microsat_data[,locus]==-9),locus] <- NA
}
missing_data <- which(is.na(microsat_data), arr.ind = T)
missing_data[,2] <- missing_data[,2]-2 # to remove the first two columns

# get number of loci and names
locus_names <- names(microsat_data)[-(1:2)]
num_of_loci <- length(locus_names)

# check size and frequency of alleles
if(check_data){
  for (locus in seq_along(locus_names)){
    print(count(microsat_data, locus_names[locus]))
  }
}

# creates vector with the position of individuals belonging to the two populations
pop1  <- which(microsat_data$population==1)
pop2  <- which(microsat_data$population==2)

# calculates sample sizes for each population (and total)
sample_size_pop1  <- length(pop1)
sample_size_pop2  <- length(pop2)
sample_size_total <- sample_size_pop1 + sample_size_pop2


cat("\n Data loaded and checked\n\n")

#########################################
# TARGET SUMMARY STATISTICS
#########################################
#if(!simulations_only){
  # Calculate summary statistics for observed data
  target_sumstats <- calculate_summary_statistics(data=microsat_data[,-c(1,2)],pop1,pop2)
  statistics_head <- colnames(target_sumstats)
  
  # write header of file
  write.table( t(statistics_head),
               file="results/target_sumstats.txt",sep=" ",
               quote=F,col.names=F,row.names=F,append=F)
  
  # write summary statistics to file
  write.table( target_sumstats,
               file      = "results/target_sumstats.txt",
               sep       = " ",
               quote     = FALSE,
               col.names = FALSE,
               row.names = FALSE,
               append    = TRUE)
#}


cat("\n Target summary statistics computed\n\n")
  
  
#########################################
# SIMULATION
#########################################

# write header of reference table 
paramaters_head <- c("model",      # 0:isolation; 1:isolation with migration
                     "theta1",     # theta1        (pop1)
                     #"thetaF",     # theta founder (pop1)
                     "theta2",     # theta1 * x2   (pop2)
                     "thetaA",     # theta1 * xA   (ancestral!=pop2)
                     #"alpha1",     # growth rate in pop1
                     #"T1",         # time demographic change (pop1)
                     #"T2",         # time demographic change (pop2)
                     "M",          # migration
                     "TS",         # time split
                     "PGSM")       # parameter for the generalised stepwise mutation model   
write.table( t(c(paramaters_head,statistics_head)),
             file      = paste0(reftable_file,".txt"),
             sep       = " ",
             quote     = F,
             col.names = F,
             row.names = F,
             append    = F)

if(parallel_sims){
  cl <- makeCluster(num_of_threads)
  registerDoParallel(cl)
  clusterEvalQ(cl,library(pegas))
  clusterEvalQ(cl,library(mmod))
  for (batch in 1:(nsim/batch_size) ){
    print(paste(batch,"of",(nsim/batch_size)))
    ref_table <- foreach(sim=seq_len(batch_size),.combine=rbind) %dopar% do_sim(sim,batch_size,theta_min,theta_max,
                                                                                M_min,M_max,T_min,T_max,
                                                                                alpha_min,alpha_max,
                                                                                PGSM_min,PGSM_max,sample_size_total,
                                                                                sample_size_pop1,sample_size_pop2,
                                                                                num_of_loci, missing_data)
    write.table( ref_table,
                 file      = paste0(reftable_file,".txt"),
                 sep       = " ",
                 quote     = F,
                 col.names = F,
                 row.names = F,
                 append    = T)
    remove(ref_table)
    gc()
    
  }
  stopCluster(cl)
  
}else{
  ref_table <- matrix(NA,nrow=nsim,ncol=length(c(paramaters_head,statistics_head)))
  for(sim in 1:nsim){
    ref_table[sim,] <- do_sim(sim,nsim,theta_min,theta_max,
                              M_min,M_max,T_min,T_max,
                              PGSM_min,PGSM_max,sample_size_total,
                              sample_size_pop1,sample_size_pop2,
                              num_of_loci, missing_data)
    
  }
  write.table( ref_table,
               file      = paste0(reftable_file,".txt"),
               sep       = " ",
               quote     = F,
               col.names = F,
               row.names = F,
               append    = T)
  
  

}
gc()

ref_table <- read.table(paste0(reftable_file,".txt"),header=T)
#dim(ref_table)
#head(ref_table)
#summary(ref_table)
save(ref_table,file=paste0(reftable_file,".RData"))

cat("\n Simulations finished\n\n")


#########################################
# ABC
#########################################
if(!simulations_only){
  # read reference table
  load(file=paste0(reftable_file,".RData"))
  
  ref_table <- ref_table[!is.na(ref_table[,"h_JostD"]),]
  ref_table <- ref_table[!is.infinite(ref_table[,"Nei_dist"]),]
  ref_table <- ref_table[!is.na(ref_table[,"v_Phi_st"]),]
  ref_table <- ref_table[!is.na(ref_table[,"g_Phi_st"]),]
  ref_table <- ref_table[!is.na(ref_table[,"m_Beta_pop1"]),]
  ref_table <- ref_table[!is.na(ref_table[,"v_Beta_pop1"]),]
  ref_table <- ref_table[!is.na(ref_table[,"v_Beta_pop2"]),]
  ref_table <- ref_table[!is.na(ref_table[,"v_Beta_total"]),]

  save(ref_table,file=paste0(reftable_file,".RData"))
  
  #sum(is.na(ref_table))
  dim(ref_table)
  

  # read target summary statistics
  target_sumstats <- read.table("results/target_sumstats.txt",header=T)
  
  # modify model column in reference table to reflect more intuitevely the two models
  ref_table[,1] <- ref_table[,1]+1 # 1 = I  - Isolation (pure divergence)
                                   # 2 = IM - Isolation with Migration (divergence with migration)
  
  # define set of summary statistics to be used (currently: all are used)
  sumstats_names <- scan(file="results/target_sumstats.txt",what=character(),nlines=1)
  #(sumstats_names==statistics_head)
  #sumstats_names <- sumstats_names[-(31:36)]
  
  # estimate posterior probabilities of both models (rejection and logistic regression methods)
  #model_selection_result <- postpr(target  = target_sumstats[sumstats_names],
  #                                 index   = ref_table[,1],
  #                                 sumstat = ref_table[,sumstats_names],
  #                                 tol     = tolerance,
  #                                 method  = "mnlogistic")
  #summary(model_selection_result)
  
  #  estimate posterior probabilities of both models (random forest)
  modindex <- as.factor(ref_table$model)
  sumsta   <- ref_table[,sumstats_names]
  model_RF <- abcrf(modindex~.,
                    data=data.frame(modindex, sumsta),
                    ntree=2000,
                    paral=T,
                    ncores=num_of_threads)
  model_RF$prior.err
  plot(model_RF,
       training=data.frame(modindex, sumsta),
       obs=target_sumstats)
  
  
  err.abcrf(model_RF,
            training=data.frame(modindex, sumsta),
            paral=T,
            ncores=num_of_threads)
  model_selection_result_RF <- predict(object         = model_RF,
                                       obs            = rbind(target_sumstats,target_sumstats),
                                       training       = data.frame(modindex, sumsta),
                                       ntree          = 1000,
                                       paral          = T,
                                       ncores         = num_of_threads,
                                       paral.predict  = T,
                                       ncores.predict = num_of_threads)
  (model_selection_result_RF)

  
  save(model_RF,
       model_selection_result_RF,
       file="results/model_choice.RData")
  remove(model_RF,
         model_selection_result_RF)
  
  #load(file="results/model_choice.RData")
  # separate referece table into two reference tables (for the two models)
  ref_tableIM <- ref_table[which(ref_table[,1]==2),]
  #remove(ref_table)
  gc()
  dim(ref_tableIM)
  
  # ISOLATION WITH MIGRATION
  #--------------------------

  #logit    <- function(x){log(x/(1-x))}
  #invlogit <- function(x){exp(x)/(1+exp(x))}

  sumstaIM   <- ref_tableIM[,sumstats_names]
  sumsta     <- ref_table[,sumstats_names]
  
  # theta1
  #------------
  
  log10theta1 <- log10(ref_table[,"theta1"])
  RFmodel_theta1 <- regAbcrf(formula = log10theta1~.,
                             data    = data.frame(log10theta1, sumsta),
                             ntree   = 1000,
                             paral   = T,
                             ncores  = num_of_threads)
  plot(x     = RFmodel_theta1,
       n.var = ncol(sumsta))
  err.regAbcrf(object   = RFmodel_theta1,
               training = data.frame(log10theta1, sumsta),
               paral    = T,
               ncores   = num_of_threads)
  posterior_theta1 <- predict(object    = RFmodel_theta1,
                              obs       = rbind(target_sumstats,target_sumstats),
                              training  = data.frame(log10theta1, sumsta),
                              quantiles = c(0.025,0.5,0.975),
                              paral     = T,
                              ncores    = num_of_threads)
  (posterior_theta1)
  densityPlot(object    = RFmodel_theta1,
              obs       = rbind(target_sumstats,target_sumstats),
              training  = data.frame(log10theta1, sumsta),
              main      = expression(log[10]*theta["W"]),
              paral     = T, 
              ncores    = num_of_threads)
  lines(density(log10theta1), col="grey")
  
  save(RFmodel_theta1,
       posterior_theta1,
       file="results/posterior_theta1.RData")
  remove(RFmodel_theta1,
         posterior_theta1)
  gc()
  
  
  
    
  
  
  # theta2
  #------------
  
  log10theta2 <- log10(ref_table[,"theta2"])
  RFmodel_theta2 <- regAbcrf(formula = log10theta2~.,
                             data    = data.frame(log10theta2, sumsta),
                             ntree   = 1000,
                             paral   = T,
                             ncores  = num_of_threads)
  plot(x     = RFmodel_theta2,
       n.var = ncol(sumsta))
  err.regAbcrf(object   = RFmodel_theta2,
               training = data.frame(log10theta2, sumsta),
               paral    = T,
               ncores   = num_of_threads)
  posterior_theta2 <- predict(object    = RFmodel_theta2,
                              obs       = rbind(target_sumstats,target_sumstats),
                              training  = data.frame(log10theta2, sumsta),
                              quantiles = c(0.025,0.5,0.975),
                              paral     = T,
                              ncores    = num_of_threads)
  (posterior_theta2)
  densityPlot(object    = RFmodel_theta2,
              obs       = rbind(target_sumstats,target_sumstats),
              training  = data.frame(log10theta2, sumsta),
              main      = expression(log[10]*theta["E"]),
              paral     = T, 
              ncores    = num_of_threads)
  lines(density(log10theta2), col="grey")
  
 
  save(RFmodel_theta2,
       posterior_theta2,
       file="results/posterior_theta2.RData")
  remove(RFmodel_theta2,
         posterior_theta2)
  gc()
  
  
  
  
  
  
   
  
  
  # thetaA
  #------------
  
  log10thetaA <- log10(ref_table[,"thetaA"])
  RFmodel_thetaA <- regAbcrf(formula = log10thetaA~.,
                             data    = data.frame(log10thetaA, sumsta),
                             ntree   = 1000,
                             paral   = T,
                             ncores  = num_of_threads)
  plot(x     = RFmodel_thetaA,
       n.var = ncol(sumsta))
  err.regAbcrf(object   = RFmodel_thetaA,
               training = data.frame(log10thetaA, sumsta),
               paral    = T,
               ncores   = num_of_threads)
  posterior_thetaA <- predict(object    = RFmodel_thetaA,
                              obs       = rbind(target_sumstats,target_sumstats),
                              training  = data.frame(log10thetaA, sumsta),
                              quantiles = c(0.025,0.5,0.975),
                              paral     = T,
                              ncores    = num_of_threads)
  (posterior_thetaA)
  densityPlot(object    = RFmodel_thetaA,
              obs       = rbind(target_sumstats,target_sumstats),
              training  = data.frame(log10thetaA, sumsta),
              main      = expression(log[10]*theta["A"]),
              paral     = T, 
              ncores    = num_of_threads)
  lines(density(log10thetaA), col="grey")
  
  
  
  save(RFmodel_thetaA,
       posterior_thetaA,
       file="results/posterior_thetaA.RData")
  remove(RFmodel_thetaA,
         posterior_thetaA)
  gc()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # TS
  #------------
  
  log10TS <- log10(ref_table[,"TS"])
  RFmodel_TS <- regAbcrf(formula = log10TS~.,
                         data    = data.frame(log10TS, sumsta),
                         ntree   = 1000,
                         paral   = T,
                         ncores  = num_of_threads)
  plot(x     = RFmodel_TS,
       n.var = ncol(sumsta))
  err.regAbcrf(object   = RFmodel_TS,
               training = data.frame(log10TS, sumsta),
               paral    = T,
               ncores   = num_of_threads)
  posterior_TS <- predict(object    = RFmodel_TS,
                          obs       = rbind(target_sumstats,target_sumstats),
                          training  = data.frame(log10TS, sumsta),
                          quantiles = c(0.025,0.5,0.975),
                          paral     = TRUE,
                          ncores    = num_of_threads)
  (posterior_TS)
  densityPlot(object    = RFmodel_TS,
              obs       = rbind(target_sumstats,target_sumstats),
              training  = data.frame(log10TS, sumsta),
              main      = expression(log[10]*"T"["F"]),
              paral     = TRUE, 
              ncores    = num_of_threads)
  lines(density(log10TS), col="grey")
  
  log10TS    <- log10(ref_tableIM[,"TS"])
  RFmodel_TS_IM <- regAbcrf(formula = log10TS~.,
                            data    = data.frame(log10TS, sumstaIM),
                            ntree   = 1000,
                            paral   = T,
                            ncores  = num_of_threads)
  plot(x     = RFmodel_TS_IM,
       n.var = ncol(sumstaIM))
  err.regAbcrf(object   = RFmodel_TS_IM,
               training = data.frame(log10TS, sumstaIM),
               paral    = T,
               ncores   = num_of_threads)
  posterior_TS_IM <- predict(object    = RFmodel_TS_IM,
                          obs       = rbind(target_sumstats,target_sumstats),
                          training  = data.frame(log10TS, sumstaIM),
                          quantiles = c(0.025,0.5,0.975),
                          paral     = TRUE,
                          ncores    = num_of_threads)
  (posterior_TS_IM)
  densityPlot(object    = RFmodel_TS_IM,
              obs       = rbind(target_sumstats,target_sumstats),
              training  = data.frame(log10TS, sumstaIM),
              main      = expression(log[10]*"T"["F"]),
              paral     = TRUE, 
              ncores    = num_of_threads)
  lines(density(log10TS), col="grey")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # M
  #------------
  
  log10M <- log10(ref_tableIM[,"M"])
  RFmodel_M <- regAbcrf(formula = log10M~.,
                        data    = data.frame(log10M, sumstaIM),
                        ntree   = 1000,
                        paral   = T,
                        ncores  = num_of_threads)
  plot(x     = RFmodel_M,
       n.var = ncol(sumstaIM))
  
  err.regAbcrf(object   = RFmodel_M,
               training = data.frame(log10M, sumstaIM),
               paral    = T,
               ncores   = num_of_threads)
  posterior_M <- predict(object    = RFmodel_M,
                         obs       = rbind(target_sumstats,target_sumstats),
                         training  = data.frame(log10M, sumstaIM),
                         quantiles = c(0.025,0.5,0.975),
                         paral     = TRUE,
                         ncores    = num_of_threads)
  (posterior_M)
  densityPlot(object    = RFmodel_M,
              obs       = rbind(target_sumstats,target_sumstats),
              training  = data.frame(log10M, sumstaIM),
              main      = expression(log[10]*"M"),
              paral     = TRUE, 
              ncores    = num_of_threads)
  lines(density(log10M), col="grey")
  

  
  
  
  
  # PGSM
  #------------
  
  PGSM <- ref_table[,"PGSM"]
  RFmodel_PGSM <- regAbcrf(formula = PGSM~.,
                           data    = data.frame(PGSM, sumsta),
                           ntree   = 1000,
                           paral   = T,
                           ncores  = num_of_threads)
  plot(x     = RFmodel_PGSM,
       n.var = ncol(sumsta))
  err.regAbcrf(object   = RFmodel_PGSM,
               training = data.frame(PGSM, sumsta),
               paral    = T,
               ncores   = num_of_threads)
  posterior_PGSM <- predict(object    = RFmodel_PGSM,
                            obs       = rbind(target_sumstats,target_sumstats),
                            training  = data.frame(PGSM, sumsta),
                            quantiles = c(0.025,0.5,0.975),
                            paral     = TRUE,
                            ncores    = num_of_threads)
  (posterior_PGSM)
  densityPlot(object    = RFmodel_PGSM,
              obs       = rbind(target_sumstats,target_sumstats),
              training  = data.frame(PGSM, sumsta),
              main      = expression("P"["GSM"]),
              paral     = TRUE, 
              ncores    = num_of_threads)
  lines(density(PGSM), col="grey")
  
  
  log10M <- log10(ref_tableIM[,"M"])
  log10TS <- log10(ref_tableIM[,"TS"])
  cov_M_TS <- covRegAbcrf(regForest1     = RFmodel_M,
                          regForest2     = RFmodel_TS_IM,
                          obs            = rbind(target_sumstats,target_sumstats),
                          training1      = data.frame(log10M, sumstaIM),
                          training2      = data.frame(log10TS, sumstaIM),
                          ntree          = 1000,
                          paral          = T,
                          ncores         = num_of_threads,
                          paral.predict  = T,
                          ncores.predict = num_of_threads)
  (cov_M_TS)
  
  
  # mutation rate (Apodemus agrarius)
  #---------------
  
  # upper limit for A. agrarius (Baker et al. 2008): 1.4E-3
  # upper limit for Apodemus sylvaticus (Baker et al. 2008): 4.7E-3
  # Mus musculus (Dietrich et al 1992) 9/(22*121*75)=4.507889e-05
  # assumed generation time = 0.5
  
  NW <- 5.25/(4 * 5E-5)      #    theat=4 N u
  NE <- 11.86/(4 * 5E-5)
  t  <- 0.17 * 4 * NW * 0.5
  m  <- 3.64 / (4 * NW)
    
    
  save(RFmodel_theta1,
       posterior_theta1,
       RFmodel_theta2,
       posterior_theta2,
       RFmodel_thetaA,
       posterior_thetaA,
       RFmodel_TS,
       posterior_TS,
       RFmodel_TS_IM,
       posterior_TS_IM,
       RFmodel_M,
       posterior_M,
       RFmodel_PGSM,
       posterior_PGSM,
       cov_M_TS,
       file="results/parameter_estimate.RData")
  #load(file="results/parameter_estimate.RData")

  # Plot Posterior probability distribution (black) compared to prior probability distribution (grey)
  
  png(filename="results/Prior_vs_Posterior.png", width=8.3, height=11.7, units = "in", res=700)
  #pdf(file="results/Prior_vs_Posterior.pdf", width=8.3, height=11.7)
  par(cex.axis = 1.5,
      cex.lab  = 2,
      cex.main = 2,
      mar      = c(2,2,3,2),
      oma      = c(2,5,0,0),
      mfrow    = c(3,2))
  
  densityPlot(object    = RFmodel_theta1,
              obs       = as.data.frame(target_sumstats),
              training  = data.frame(log10theta1, sumsta),
              main      = expression(log[10]*italic(theta)["W"]),
              paral     = T, 
              ncores    = num_of_threads)
  lines(density(log10theta1), col="grey")
  
  densityPlot(object    = RFmodel_theta2,
              obs       = as.data.frame(target_sumstats),
              training  = data.frame(log10theta2, sumsta),
              main      = expression(log[10]*italic(theta)["E"]),
              paral     = T, 
              ncores    = num_of_threads)
  lines(density(log10theta2), col="grey")
  
  densityPlot(object    = RFmodel_thetaA,
              obs       = as.data.frame(target_sumstats),
              training  = data.frame(log10thetaA, sumsta),
              main      = expression(log[10]*italic(theta)["A"]),
              paral     = T, 
              ncores    = num_of_threads)
  lines(density(log10thetaA), col="grey")

    densityPlot(object    = RFmodel_TS,
              obs       = as.data.frame(target_sumstats),
              training  = data.frame(log10TS, sumsta),
              main      = expression(log[10]*italic("T")["F"]),
              paral     = TRUE, 
              ncores    = num_of_threads)
  lines(density(log10TS), col="grey")
  
  densityPlot(object    = RFmodel_M,
              obs       = as.data.frame(target_sumstats),
              training  = data.frame(log10M, sumstaIM),
              main      = expression(log[10]*italic("M")),
              paral     = TRUE, 
              ncores    = num_of_threads)
  lines(density(log10M), col="grey")
  
  densityPlot(object    = RFmodel_PGSM,
              obs       = as.data.frame(target_sumstats),
              training  = data.frame(PGSM, sumsta),
              main      = expression(italic("P")["GSM"]),
              paral     = TRUE, 
              ncores    = num_of_threads)
  lines(density(PGSM), col="grey")
  
  plot.new()
  legend(x      = "topleft",
         legend = c("Prior probability distribution", "Posterior probability distribution"),
         lty    = 1,
         cex    = 1.5,
         col    = c("grey","black"))
  mtext(text  = "Probability density",
        side  = 2,
        line  = 1,
        adj   = 0.5,
        cex   = 1.5,
        outer = TRUE)
  
  dev.off ( which=dev.cur() )
  
  
  
  
    
}  


