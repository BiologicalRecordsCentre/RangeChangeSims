# Nick Isaac
# with Arco van Strien
# Initial code written 11-16 July 2012
# Simulations to compare the performance of range change methods
#
###################################################################################################################
#
# TO DO
# Extract mean trend estimates (esp under power test): Frescalo probably under-estimates slope
#
###################################################################################################################
#
# ARCO & MARNIX - Running this script should work in Windows so long as:
# 1. you define your local Working directory (delete the code down to the next refenrece to your name (line 52 in this version)
# 2. the file 'Sim function 121022.r' is stored in that working directory
# 3. if saving copies of the raw data (variable sv, below) then you need a subfolder called 'SimRecords'
#
###################################################################################################################
# VERSION HISTORY
# 24-25 April: I changed the way sites are selected to be visited to the power-law (nVisits is no longer a parameter).
#   Added the Maes et al method.
#   simdata files (the raw records) now contain information about the parameters in the name.
# 8-9 April: Various changes agreed during meeting in Netherlands 2-5 April. Principle ones are:
#   Disaggregated iterate_all_scenarios() into two functions
#   simpars is now an attribute of records (as well as output)
#   Removed the tendency for speciose sites to be visited more often
#   Replaced pF with separate Occ and DetP (in a list)
#   Replaced the 'visit-based' mixed model with the 'Binomial' version (same formulation but quicker and more stable)
# 7 -13 Feb: added greater user control of parameters
# 25 Jan - 4 Feb: Added B3 & F scenarios
# 7-17 January
#   Various modification of input parameters to explore performance of Frescalo
# 30 November - 3 December: Reparameterization
#   create_data uses a beta distribution for occupancy and also includes detection probabilities
#   species detection probabilities now used in 
#   removed some scenarios
#   altered detail of incomplete recording (balance, starting values and final)
#   changed how sites are selected, so a greater proportion receive multiple visits
# 20 November: Added Frescalo (NB only works with 1000 sites)
#   changed test of power to 30% decline (IUCN threshold for VU under A2c) from 25%
#   removed some of the combination scenarios
#   added a new scenario E: large numbers of incidental records of just one species
#   fixed error in estimation of power
## 3 November: I modified iterate_all_scenarios to allow multiple thresholds for the Mixed model
#	Also I added a resample() function into shorten_lists() to avoid the problem of sampling from a vector of length 1.
## 25 October: I started completely changing the sim structure. It's much quicker this way
# 	Now I run a single recording scenario then subset the data according to rules in order to generate the alternatives.
#	Also, fitting the models has been streamlined and quickened
#	Added a test of power
## 16-17 October: copied. Fixed a few errors in the functions. Otherwise dabbled (changing num_reps here)
# 13 July: runs 7 scenarios of recording behaviour with 6 methods for estimating change
#
# This version is compatible with Linux
#
###################################################################################################################
###################################################################################################################
###################################################################################################################
rm(list=objects())
#options(warn=2) # sets warnings to errors

#set the working directory
if (grepl("linux", R.version$platform)){
    if(grepl("redhat", R.version$platform)){    #WLLFXXX
        homedir<-'/prj/NEC04273/Sims' # 'prj' changed from '~' 17 Jan
    } else {#neodarwin
        homedir <- '~/myR'
    }
} else homedir <- 'P:/NEC04273_SpeciesDistribution/Workfiles/Range change sims' # same location as 'Sim functions 12xxxx.r'
setwd(homedir)

#point the machine to where the libraries are stores
if(any(grepl("redhat", R.version$platform), grepl("w32", R.version$platform))){    #WLLF006 or Windoze
    libloc <- NULL
} else libloc <- 'packages' # neodarwin

library(lme4, lib=libloc) # this is not necessary on Windows - the packages are called in specific functions 
library(reshape2, lib=libloc) # this is not necessary on Windows - the packages are called in specific functions

###################################################
#
# ARCO & MARNIX - replace all code above with simple reference to the local Working Directory
# 
################################################################################################################### 
datecode <- format(Sys.Date(),'%y%m%d')
source('Sim_functions.r')
################################################################################################################### FULL SIMS
# Set the parameters here, to be passed further down into the workhorse functions
nSi=1000; nY=10
ps=list(init=0.6, final=0.9) # from c(1/3, 2/3): higher starting value but less steep increase
di=0.2 # parameters defining the degree of bias and under-recording
sv = F # save data files - set to TRUE to save copies of the 'records' (adds very little time)
pF = list(Occ=0.5, DetP=0.5) # the probability of occurrence & detection for the focal species
d= 0 # extinction rate of the focal species (over the full timeframe, not per annum)
sd=F # should the data be saved to
nSp=25
pSVS=0.05 # the proportion of sites receiving a single visit
vrs=list(sel=FALSE,num=TRUE) # correlation between richness and probability of being a) selected for a visited and b) number of visits
mv = 10 #maximum number of visits to any one site in any one year
st=T # should the number of visits vary stochastically from year to year? Say yes, since this is what i used before (inadvertently)
combos=F
Sc <- 'BCDF' # which scenarios should be run. At least one of B-F must be selected
Fr = 2#:1 # should the frescalo method be run - set to FALSE to bypass this
MM = c(2,4)#2:4 # should the Mixed model be run - set to FALSE to bypass this (about 50% of the total runtime excl Frescalo)
nr = 500

################################################################################################################### EXPLORE
#true_data <- create_data(nSites=1000, nSpecies=25)
#records <- generate_records(nYrs=10, pSVS=pSVS, true_data) # multiple years
records <- generate_records(nYrs=10, pSVS=pSVS, true_data=create_data()) # multiple years
#records <- generate_all_scenarios(nSites=nSi, nSpecies=nSp, nYrs=nY, pSVS=pSVS, pFocal=pF, p_short=ps, pDetMod=di, decline=d, Scenarios=Sc)
#output <- iterate_all_scenarios(nreps=nr, nSpecies=nSp, nSites=nSi, nYrs=nY, pSVS=pSVS, p_short=ps, pDetMod=di, 
#                                pFocal=pF, decline=d, id=code, save_data=sv, inclMM=MM, Frescalo=Fr, Scenarios=Sc)
################################################################################################################### EXPLORE

################
# MARNIX - just change this next line to TRuE
ExtractOnly = F
################

### temp

#Sc='D' # one other scenario - the output querying doesn't work with only A
Fr=0
MM=0
nSi=1000
#nr=20#50
d=0.3
#pSVS=0.05 # the proportion of sites receiving a single visit
#nr=1
system.time({
    for(d in c(0, 0.3)) 
    for(pSVS in c(0.05, 0.1, 0.07)) 
        {
        # prepare the run
        ch <- ifelse(d==0, 'V', 'P')
        code = paste(ch,'_',pSVS*100,'SVS_',datecode, sep='') 
        
        if(nSi == 1000) set.seed(153) #for continuity, but only for comparing Dutch & British datasets 
        else Fr <- F
        
        if (!ExtractOnly) {
            # wE ARE IN THE UK, running the trend methods
            output <- iterate_all_scenarios(nreps=nr, nSpecies=nSp, nSites=nSi, nYrs=nY, pSVS=pSVS, p_short=ps, pDetMod=di, 
                                            mv=mv, vrs=vrs, stoch=st, Scenarios=Sc, combos=combos, pFocal=pF, decline=d, id=code, save_data=sd,
                                            inclMM=MM, Frescalo=Fr)
            get_all_stats(output, save_to_txt=T) #single output including error rates and stats
       } else { 
            # We are in the Netherlands running Occupancy
           replicate(nr, recs <- generate_all_scenarios(nSpecies=nSp, nSites=nSi, nYrs=nY, pSVS=pSVS, p_short=ps, pDetMod=di, 
                                          mv=mv, vrs=vrs, stoch=st,Scenarios=Sc, pFocal=pF, combos=combos, decline=d, id=code, save_data=T))
       }
    }
})
# end

