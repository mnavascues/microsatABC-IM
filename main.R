# Script to perform an ABC analysis
# Miguel Navascués 2016
# Miguel.Navascues@supagro.inra.fr

#########################################
# ABC GENERAL SETTINGS
#########################################

nsim          <- 100000                   # number of simulations
reftable_file <- "reference_table.txt"    # reference table file name
tolerance     <- 0.01                     # tolerance level (proportion of kept simulations)

# Model description
#
#    pop1              pop2
#
#   theta1            theta2
#  \      /    M       /  \
#   \    /    <--     /    \
#    \  /     -->    /      \
#     \/            /        \
#  TS ------------------------
#          |     |
#          |     |
#          |     |
#          |     |
#          thetaA
#


#########################################
# DEFINE PRIORS (upper and lower limit)
#########################################

theta_min <- -3
theta_max <-  3

alpha_min  <- 0
alpha_max  <- 10

M_min <- 0
M_max <- 100

T_min <- 0
T_max <- 10

PGSM_min <- 0     # Generalised Stepwise Mutation model parameter 
PGSM_max <- 1     # (probability of mutations larger than a single repeat unit)

#########################################
# LOAD REQUIRED FUNCTIONS/PACKAGES
#########################################

# load required libraries
library(plyr) # for function count() that helps to explore allele size distribution
library(pegas)  # for calculating heterozygosity: function H()
library(abc)    # for approximate Bayesian computation functions
require(hexbin) # for plots of PCA
require(grid)   # for plots of PCA
require(laeken) # necessary? check if it is used...
require(abcrf)  # abc random forest

# Uses RR Hudson coalescent simulator (ms)
# load RR Hudson's R script distributed with ms to read ms output in R
# http://home.uchicago.edu/~rhudson1/source/mksamples.html
source("readms.output.R")

# load other functions
source("fun.R")

#########################################
# LOAD DATA
#########################################

# read real data file
microsat_data <- read.table("agrarius.str",header=T)

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
for (locus in seq_along(locus_names)){
  print(count(microsat_data, locus_names[locus]))
}

# creates vector with the position of individuals belonging to the three populations
pop1  <- which(microsat_data$POP==1)
pop2  <- which(microsat_data$POP==2)

# calculates sample sizes for each population (and total)
sample_size_pop1  <- length(pop1)
sample_size_pop2  <- length(pop2)
sample_size_total <- sample_size_pop1 + sample_size_pop2
num_of_loci <- ncol(microsat_data)-2

#########################################
# TARGET SUMMARY STATISTICS
#########################################

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

#########################################
# SIMULATION
#########################################

# write header of reference table 
paramaters_head <- c("model",      # 0:isolation; 1:isolation with migration
                     "theta1",     # theta1        (pop1)
                     "theta2",     # theta1 * x2   (pop2)
                     "thetaA",     # theta1 * xA   (ancestral)
                     "alpha1",     # growth rate   (pop1)
                     "alpha2",     # growth rate   (pop1)
                     "mig12",      # migration
                     "mig21",      # migration
                     "TS",         # time split
                     "PGSM")       # parameter for the generalised stepwise mutation model   
write.table( t(c(paramaters_head,statistics_head)),
             file      = reftable_file,
             sep       = " ",
             quote     = F,
             col.names = F,
             row.names = F,
             append    = F)

for(sim in 1:nsim){
  # write progress of simulation on screen 
  if (sim==1)      cat(paste(sim,"of",nsim))
  if (sim%%100==0) cat(paste(sim,"of",nsim))
  if (sim==nsim)   cat(paste(sim,"of",nsim))
  
  # SAMPLE FROM PRIOR
  #-------------------
  
  # model (migration OR pure divergence) is sampled from prior
  model <- sample(c(0,1),1)
  
  # take parameter values from priors
  theta1 <- 10^runif( 1, min=theta_min, max=theta_max )
  theta2 <- 10^runif( 1, min=theta_min, max=theta_max )
  thetaA <- 10^runif( 1, min=theta_min, max=theta_max )
  x2     <- theta2/theta1
  xA     <- thetaA/theta1
  alpha1 <- runif( 1, min=alpha_min, max=alpha_max)
  alpha2 <- runif( 1, min=alpha_min, max=alpha_max)
  if (model==1){
    mig12   <- runif( 1, min=M_min, max=M_max )
    mig21   <- runif( 1, min=M_min, max=M_max )
  }else{
    mig21<-mig12<-0
  }
  TS <- runif( 1, min=T_min, max=T_max )
  PGSM <- runif( 1, min=PGSM_min, max=PGSM_max )
  
  # RUM MS
  #--------
  
  # generate text with ms command
  ms_out_file <- "msout.txt"
  ms_run <- paste( "ms", sample_size_total, num_of_loci,      # total sample size & number of loci
                   "-t", theta1,                              # population size in pop1
                   "-I 2", sample_size_pop1, sample_size_pop2,# sample sizes per population
                   "-n 2", x2,                                # population size in pop2
                   "-m 1 2", mig12,                 # migration rate from pop1 to pop2
                   "-m 2 1", mig21,
                   "-g 1", alpha1,                   # growth rate in pop1
                   "-g 2", alpha2,                   # growth rate in pop1
                   "-ej", TS, "1 2",                # creation of pop1 from pop2
                   "-en", TS, "2", xA,              # population size in ancestral
                   ">", ms_out_file)                # output file
  
  # rum ms on system
  if(.Platform$OS.type == "unix") {
    ms_run <- paste( "./", ms_run, sep="" )
    system( ms_run )
  }
  
  # READ SIMULATED DATA AND TRANSFORM TO MICROSATELLITE
  #-----------------------------------------------------
  
  # read simulation output (generated with an infinite site model, ISM)
  msout <- read.ms.output(file=ms_out_file)
  
  # transform simulated data from ISM into microsatellite evolving under GSM
  # (generalised stepwise mutation model)
  mutations <- msout$segsites
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

  sim_sumstats <- calculate_summary_statistics(data=simulated_data,pop1,pop2)

  # WRITE PARAMATERS AND SUMMARY STATISTICS TO REFERENCE TABLE
  #------------------------------------------------------------
  write.table( cbind( ## PARAMETERS
                      model,      # 0:isolation; 1:isolation with migration
                      theta1,     # theta1 (pop1)
                      theta2,     # theta2 (pop2)
                      thetaA,     # thetaA (ancestral)
                      alpha1,      # growth rate (pop1)
                      alpha2,      # growth rate (pop2)
                      mig12,      # migration
                      mig21,      # migration
                      TS,         # time split of pop1 from pop2
                      PGSM,
                      #### SUMMARY STATISTICS                      
                      sim_sumstats),
               file=reftable_file,sep=" ",
               quote=F,col.names=F,row.names=F,append=T)
  
} # end for(sim in 1:nsim)

#########################################
# ABC: REJECTION + REGRESSION
#########################################

# read reference table
ref_table <- read.table("reference_table.txt",header=T)
#dim(ref_table)
#head(ref_table)
#summary(ref_table)

# read target summary statistics
target_sumstats<-read.table("target_sumstats.txt",header=T)

# modify model column in reference table to reflect more intuitevely the two models
ref_table[which(ref_table[,1]==0),1] <- "I"    # Isolation (pure divergence)
ref_table[which(ref_table[,1]==1),1] <- "IM"   # Isolation with Migration (divergence with migration)

# define set of summary statistics to be used (currently: all are used)
sumstats_names <- scan(file="target_sumstats.txt",what=character(),nlines=1)
(sumstats_names==statistics_head)
sumstats_names <- sumstats_names[-(31:36)]

# estimate posterior probabilities of both models (rejection and logistic regression methods)
model_selection_result <- postpr(target  = target_sumstats[sumstats_names],
                                 index   = ref_table[,1],
                                 sumstat = ref_table[,sumstats_names],
                                 tol     = tolerance,
                                 method  = "mnlogistic")
summary(model_selection_result)

#  estimate posterior probabilities of both models (random forest)
for( no_sims in c(10000,20000,30000,40000,50000,60000) ){
  model_RF <- abcrf(modindex = as.factor(ref_table$model)[seq_len(no_sims)],
                          sumsta   = ref_table[seq_len(no_sims),sumstats_names])
  assign(paste("model_RF",no_sims,sep="_"),model_RF )
  remove(model_RF)
}
model_RF_10000$prior.err
model_RF_20000$prior.err
model_RF_30000$prior.err
model_RF_40000$prior.err
model_RF_50000$prior.err
model_RF_60000$prior.err

err.abcrf(model_RF_40000)

model_selection_result_RF <- predict(object   = model_RF_40000,
                                     obs      = target_sumstats[sumstats_names],
                                     sampsize = 40000,
                                     ntree    = 500,
                                     paral    = TRUE)
(model_selection_result_RF)
summary(model_selection_result_RF)

# separate referece table into two reference tables (for the two models)
ref_tableIM <- ref_table[which(ref_table[,1]=="IM"),]
ref_tableI  <- ref_table[which(ref_table[,1]=="I") ,]
remove(ref_table)

# ISOLATION WITH MIGRATION
#--------------------------

# estimate posterior probability distribution of parameter values
abcresultIM <- abc(target  = target_sumstats[sumstats_names],
                   param   = ref_tableIM[,c(2:6,8:9)],
                   sumstat = ref_tableIM[,sumstats_names],
                   tol     = tolerance,
                   method  = "loclinear",
                   transf  = c("log","log","log",
                               "none","none","none",
                               "logit"),
                   logit.bounds=matrix(c(PGSM_min,PGSM_max),nrow=7,ncol=2,byrow=T))
summary(abcresultIM)

# calculate principal components for summary statistics from simulations
PCA_stats  <- princomp(ref_tableIM[,sumstats_names])
# predict location of real data in PC space
PCA_target <- predict(PCA_stats, target_sumstats[sumstats_names])
# represent graphically PCA to verify that simulations are producing 
# simulations similar to real data
pdf(file="PCA_IM.pdf", width=11.7, height=8.3)
for (pci in 1:7){
  for (pcj in (pci+1):8){
    hbin<-hexbin(PCA_stats$scores[,pci], PCA_stats$scores[,pcj],xbins=100,xlab=paste("PC",pci),ylab=paste("PC",pcj))
    pp<-plot(hbin,legend=FALSE, main="PCA on summary statistics")
    pushHexport(pp$plot.vp)
    grid.points(PCA_target[pci],PCA_target[pcj],pch="*",gp=gpar(cex=2,col="red"))
  }
}
dev.off ( which=dev.cur() )

# ISOLATION WITH MIGRATION
#--------------------------

# estimate posterior probability distribution of parameter values
abcresultI <- abc(target  = target_sumstats[sumstats_names],
                  param   = ref_tableI[,c(2:5,8:9)],
                  sumstat = ref_tableI[,sumstats_names],
                  tol     = tolerance,
                  method  = "loclinear",
                  transf  = c("log","log","log",
                              "none","log",
                              "logit"),
                  logit.bounds=matrix(c(PGSM_min,PGSM_max),nrow=6,ncol=2,byrow=T))
summary(abcresultI)


# calculate principal components for summary statistics from simulations
PCA_stats  <- princomp(ref_tableI[no_NA,sumstats_names])
# predict location of real data in PC space
PCA_target <- predict(PCA_stats, target_sumstats[sumstats_names])
# represent graphically PCA to verify that simulations are producing 
# simulations similar to real data
pdf(file="PCA_I.pdf", width=11.7, height=8.3)
for (pci in 1:7){
  for (pcj in (pci+1):8){
    hbin<-hexbin(PCA_stats$scores[,pci], PCA_stats$scores[,pcj],xbins=100,xlab=paste("PC",pci),ylab=paste("PC",pcj))
    pp<-plot(hbin,legend=FALSE, main="PCA on summary statistics")
    pushHexport(pp$plot.vp)
    grid.points(PCA_target[pci],PCA_target[pcj],pch="*",gp=gpar(cex=2,col="red"))
  }
}
dev.off ( which=dev.cur() )


