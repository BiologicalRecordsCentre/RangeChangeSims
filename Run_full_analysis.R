# Script to reproduce results in paper 
# Tom August
# 03/10/2013

# WARNING: This script will install libraries needed for the analysis.
# This includes updateing lme4 if the version you have is <1.0. If you
# do not want to have these packages installed do not run this script.

# install required packages if not installed
req_pkgs <- c('lme4','reshape2','abind','sp','gdata','lattice','Matrix')
inst_pkgs <- req_pkgs[!req_pkgs %in% installed.packages()]
if(length(inst_pkgs) > 0){
  install.packages(inst_pkgs)
}

if(!'sparta' %in% installed.packages()){
  install.packages('sparta_0.1.19.zip', repos = NULL)
}

# Check that we have a version of lme4 > 1.0
lme4_details <- installed.packages()['lme4',]
if(!lme4_details['Version']>1.0){
  install.packages('lme4')
}

library(sparta)

# Create the data
source('Sim_Wrapper.r')

