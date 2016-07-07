#!/usr/bin/env Rscript

# Script to perform an ABC analysis
# Miguel Navascués 2016
# Miguel.Navascues@supagro.inra.fr

### Remove last features
rm(list=ls())
ls()


#########################################
# ABC GENERAL SETTINGS
#########################################

nsim             <- 1000000             # number of simulations
input_filename   <- "agrarius.str"      # data file
reftable_file    <- "reference_table"   # reference table file name
tolerance        <- 0.01                # tolerance level (proportion of kept simulations)
seed             <- 1234 ; 
set.seed(seed,"Mersenne-Twister")
check_data       <- T
simulations_only <- F
parallel_sims    <- T
num_of_threads   <- 8

# Model graphical description
#
#     pop1            pop2
#
#    theta1          theta2
#  |       |    M     | | 
#  |__   __|  <--->   | | 
#     | |            _| |_
#     | |           |     |
#  TS |_|___________|     |
#                   |     |
#                   |     |
#                   |     |
#                   |     |
#                   thetaA
#
# NB: a population expansion is assumed in pop 1,
# but no assumption on thhe type of demographic change
# is made for pop 2 (see file fun.R to check ms command)

#########################################
# DEFINE PRIORS (upper and lower limit)
#########################################

theta_min <- -1
theta_max <-  3

M_min <- -3
M_max <- 2

T_min <- 0
T_max <- 10

PGSM_min <- 0     # Generalised Stepwise Mutation model parameter 
PGSM_max <- 1     # (probability of mutations larger than a single repeat unit)

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
source("readms.output.R")

# load other functions (from file distributed together with present file)
source("fun.R")

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
num_of_loci       <- ncol(microsat_data)-2


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
               file="target_sumstats.txt",sep=" ",
               quote=F,col.names=F,row.names=F,append=F)
  
  # write summary statistics to file
  write.table( target_sumstats,
               file      = "target_sumstats.txt",
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
                     "thetaF",     # theta founder (pop1)
                     "theta2",     # theta1 * x2   (pop2)
                     "thetaA",     # theta1 * xA   (ancestral=pop2)
                     "T1",         # time demographic change (pop1)
                     "T2",         # time demographic change (pop1)
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
  ref_table <- foreach(sim=seq_len(nsim),.combine=rbind) %dopar% do_sim(sim,nsim,theta_min,theta_max,
                                                                        M_min,M_max,T_min,T_max,
                                                                        PGSM_min,PGSM_max,sample_size_total,
                                                                        sample_size_pop1,sample_size_pop2,
                                                                        num_of_loci, missing_data)
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
}
write.table( ref_table,
             file      = paste0(reftable_file,".txt"),
             sep       = " ",
             quote     = F,
             col.names = F,
             row.names = F,
             append    = T)

ref_table <- read.table(paste0(reftable_file,".txt"),header=T)
#dim(ref_table)
#head(ref_table)
#summary(ref_table)
save(ref_table,file=paste0(reftable_file,".RData"))

cat("\n Simulations finished\n\n")


#########################################
# ABC: REJECTION + REGRESSION
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
  target_sumstats <- read.table("target_sumstats.txt",header=T)
  
  # modify model column in reference table to reflect more intuitevely the two models
  ref_table[which(ref_table[,1]==0),1] <- "I"    # Isolation (pure divergence)
  ref_table[which(ref_table[,1]==1),1] <- "IM"   # Isolation with Migration (divergence with migration)
  
  # define set of summary statistics to be used (currently: all are used)
  sumstats_names <- scan(file="target_sumstats.txt",what=character(),nlines=1)
  (sumstats_names==statistics_head)
  #sumstats_names <- sumstats_names[-(31:36)]
  
  # estimate posterior probabilities of both models (rejection and logistic regression methods)
  #model_selection_result <- postpr(target  = target_sumstats[sumstats_names],
  #                                 index   = ref_table[,1],
  #                                 sumstat = ref_table[,sumstats_names],
  #                                 tol     = tolerance,
  #                                 method  = "mnlogistic")
  #summary(model_selection_result)
  
  #  estimate posterior probabilities of both models (random forest)
  for( no_sims in c(5000,10000,15000,20000,25000,30000) ){
    model_RF <- abcrf(modindex = as.factor(ref_table$model)[seq_len(no_sims)],
                      sumsta   = ref_table[seq_len(no_sims),sumstats_names])
    assign(paste("model_RF",no_sims,sep="_"),model_RF )
    remove(model_RF)
  }
  model_RF_5000$prior.err
  model_RF_10000$prior.err
  model_RF_15000$prior.err
  model_RF_20000$prior.err
  model_RF_25000$prior.err
  model_RF_30000$prior.err
  
  err.abcrf(model_RF_10000)
  err.abcrf(model_RF_25000)
  err.abcrf(model_RF_30000)
  
  plot(model_RF_10000)
  plot(model_RF_25000)
  plot(model_RF_30000)
  
  model_selection_result_RF <- predict(object   = model_RF_10000,
                                       obs      = target_sumstats[,sumstats_names],
                                       sampsize = 10000,
                                       ntree    = 500,
                                       paral    = TRUE)
  model_selection_result_RF <- predict(object   = model_RF_25000,
                                       obs      = target_sumstats[,sumstats_names],
                                       sampsize = 25000,
                                       ntree    = 500,
                                       paral    = TRUE)
  model_selection_result_RF <- predict(object   = model_RF_30000,
                                       obs      = target_sumstats[,sumstats_names],
                                       sampsize = 25000,
                                       ntree    = 500,
                                       paral    = TRUE)
  (model_selection_result_RF)
  summary(model_selection_result_RF)
  
  model_RF_25000$model.lda
  model_RF_10000$model.lda
  
  save(model_RF_5000,  model_RF_10000, model_RF_15000,
       model_RF_20000, model_RF_25000, model_RF_30000,
       model_selection_result_RF, file="model_choice.RData")
  
  load(file="model_choice.RData")
  # separate referece table into two reference tables (for the two models)
  ref_tableIM <- ref_table[which(ref_table[,1]=="IM"),]
  ref_tableI  <- ref_table[which(ref_table[,1]=="I") ,]
  remove(ref_table)
  remove(ref_tableI)
  gc()
  dim(ref_tableIM)
  
  # ISOLATION WITH MIGRATION
  #--------------------------

  num_of_trees <- 1000
  quant <- seq(0,1,0.001)
  num_of_sim <- 30000

  #logit    <- function(x){log(x/(1-x))}
  #invlogit <- function(x){exp(x)/(1+exp(x))}

  for (param in 2:length(paramaters_head)){
    paramater_values <- ref_tableIM[seq_len(num_of_sim),param]
    if (param==2 | param==3 | param == 4 | param == 5){
      paramater_values <-log10(paramater_values)
    }
    #if (param==10){
    #  paramater_values <-logit( paramater_values )
    #}
    
    RFmodel <- quantregForest(x = ref_tableIM[seq_len(num_of_sim),sumstats_names],
                              y = paramater_values,
                              ntree = num_of_trees)

    pdf(file=paste0(paramaters_head[param],".pdf"))

    plot(RFmodel,main=paramaters_head[param])
    plot(round(num_of_trees/2) : num_of_trees,
         RFmodel$mse[round(num_of_trees/2) : num_of_trees],
         type="l",
         ylab="Error",xlab="trees",
         main=paramaters_head[param])
    RFmodel$mse[num_of_trees]
    
    plot(paramater_values,RFmodel$predicted,
         ylab="predicted",xlab="true value",
         main=paramaters_head[param])

    if (param==2 | param==3 | param == 4 | param == 5){
      plot(10^(paramater_values),10^(RFmodel$predicted),
           ylab="predicted",xlab="true value",
           main=paste0(paramaters_head[param], ", natural scale"))
    }
    #if (param==10){
    #  plot(invlogit(paramater_values),invlogit(RFmodel$predicted),
    #       ylab="predicted",xlab="true value",
    #       main=paste0(paramaters_head[param], ", natural scale"))
    #}
    
    posterior <- predict(object  = RFmodel,
                         newdata = target_sumstats[,sumstats_names],
                         what    = quant )
    posterior_mean <- predict(object  = RFmodel,
                              newdata = target_sumstats[,sumstats_names],
                              what    = mean )
    
    prior      <- quantile(paramater_values, probs=quant)
    prior_mean <- mean(paramater_values)
    
    plot(prior,
         posterior,
         xlab=expression("prior quantiles"),
         ylab=expression("posterior quantiles"),
         type="l",
         main=paramaters_head[param])
    abline(h=posterior[which(quant==0.5)],col="red",lwd=2)
    abline(h=posterior_mean,col="blue",lwd=2)
    abline(h=posterior[which(quant==0.025)],col="red",lwd=2,lty=2)
    abline(h=posterior[which(quant==0.975)],col="red",lwd=2,lty=2)
    abline(h=posterior[which(quant==0.25)],col="red",lty=2)
    abline(h=posterior[which(quant==0.75)],col="red",lty=2)
    
    varImpPlot(RFmodel,main=paste("Variable Importance to estimate ",paramaters_head[param]))
    
    dev.off()
    
    
    
    
    model_name <- paste("RFreg",paramaters_head[param],sep="_")
    assign(  model_name, RFmodel)
    
    save(get(model_name),file=paste(model_name,"parameter_estimate.RData",sep="_"))
  }

  
  paramater_values <- ref_tableIM[seq_len(num_of_sim),3]/ref_tableIM[seq_len(num_of_sim),2]
  RFreg_X1 <- quantregForest(x = ref_tableIM[seq_len(num_of_sim),sumstats_names],
                               y = paramater_values,
                               ntree = num_of_trees)
  pdf(file="X1.pdf")
  plot(RFreg_X1,main="X1")
  plot(round(num_of_trees/2) : num_of_trees,
       RFreg_X1$mse[round(num_of_trees/2) : num_of_trees],
       type="l",
       ylab="Error",xlab="trees",
       main="X1")
  RFreg_X1$mse[num_of_trees]
    
  plot(paramater_values,RFreg_X1$predicted,
       ylab="predicted",xlab="true value",
       main="X1")
  posterior <- predict(object  = RFreg_X1,
                       newdata = target_sumstats[,sumstats_names],
                       what    = quant )
  posterior_mean <- predict(object  = RFreg_X1,
                            newdata = target_sumstats[,sumstats_names],
                            what    = mean )
    
  prior      <- quantile(paramater_values, probs=quant)
  prior_mean <- mean(paramater_values)
    
  plot(prior,
       posterior,
       xlab=expression("prior quantiles"),
       ylab=expression("posterior quantiles"),
       type="l",
       main="X1")
  abline(h=posterior[which(quant==0.5)],col="red",lwd=2)
  abline(h=posterior_mean,col="blue",lwd=2)
  abline(h=posterior[which(quant==0.025)],col="red",lwd=2,lty=2)
  abline(h=posterior[which(quant==0.975)],col="red",lwd=2,lty=2)
  abline(h=posterior[which(quant==0.25)],col="red",lty=2)
  abline(h=posterior[which(quant==0.75)],col="red",lty=2)
    
  varImpPlot(RFreg_X1,main=paste("Variable Importance to estimate ","X1"))
    
  dev.off()

  save(RFreg_X1,file="X1_parameter_estimate.RData")
  
  #load(file="parameter_estimate.RData")

  pdf(file="SS_obs_vs_sim_IM.pdf", width=11.7, height=8.3)
  for(ss in seq_along(sumstats_names)){
    ss_name <- sumstats_names[ss] 
    #hist(ref_tableIM[,ss_name])
    #abline(v=target_sumstats[,ss_name],col="red")
    plot(density (ref_tableIM[,ss_name]),main=ss_name)
    abline(v=target_sumstats[,ss_name],col="red")
  }
  dev.off ( which=dev.cur() )
  
}  


