# Script to reproduce results in paper 
# Tom August
# 19/02/2014

# WARNING: This script will install libraries needed for the analysis.
# This includes updateing lme4 if the version you have is <1.0. If you
# do not want to have these packages installed do not run this script.
rm(list=ls())

# install required packages if not installed
req_pkgs <- c('lme4','reshape2','abind','sp','gdata',
              'lattice','Matrix','ggplot2','Rcpp','R2jags')
inst_pkgs <- req_pkgs[!req_pkgs %in% installed.packages()]
if(length(inst_pkgs) > 0){
  install.packages(inst_pkgs)
}

# install the sparta package (included in the github files), needed for the Frescalo
# analyses
if(!'sparta' %in% installed.packages()){
  install.packages('sparta_0.1.20.zip', repos = NULL)
}
library(sparta)

# Check that we have a version of lme4 > 1.0
# If not, update
lme4_details <- installed.packages()['lme4',]
if(!lme4_details['Version']>1.0){
  install.packages('lme4')
}

# Test JAGS is functional
source('jagsTest.r')

# Run the analysis 
# The original analysis uses 1000 runs
source('Sim_Wrapper.r')
Sim_Wrapper(number_of_runs=1)

# Plot figures
plots <- Explore_results(readPath='Output',writePath='Results')