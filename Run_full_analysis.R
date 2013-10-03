# Script to reproduce results in paper 
# Tom August
# 03/10/2013

# sparta package is required
if(!'sparta' %in% installed.packages()){
  install.packages('sparta_0.1.19.zip', repos = NULL)
}

library(sparta)

# Create the data
source('Sim_Wrapper.r')

