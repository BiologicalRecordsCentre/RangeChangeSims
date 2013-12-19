# Nick Isaac
# with Arco van Strien
# Initial code written 11-16 July 2012
# Simulations to compare the performance of range change methods
Sim_Wrapper <- function(number_of_runs=NULL){
  
  app_tmin  <- (number_of_runs * 1500)/60
  app_hours <- app_tmin%/%60
  app_min <- round(app_tmin%%60)
  cat(paste('Undertaking', number_of_runs, 'runs will take about',
            app_hours, 'hours', app_min, 'minutes.\n',
            'Should be done around',format(Sys.time()+(app_tmin*60),'%H:%M'),'\n',
            sep=' '))
  if(app_tmin>30) cat("It's probably a good time to go grab a coffee")
  
  library(lme4) 
  library(reshape2)
  library(abind)
  
  datecode <- format(Sys.Date(),'%y%m%d')
  source('Sim_functions.r')
  source('sims_to_frescalo.r')
  source('run_fresc_param_sims.r')
  
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
Fr = 1:2# should the frescalo method be run - set to FALSE to bypass this
MM = c(2,4)# should the Mixed model be run - set to FALSE to bypass this (about 50% of the total runtime excl Frescalo)
nr = number_of_runs # Number of runs

  ################################################################################################################### EXPLORE
  records <- generate_all_scenarios(nSites=nSi, nSpecies=nSp, nYrs=nY, pSVS=pSVS, pFocal=pF, p_short=ps, pDetMod=di, decline=d, Scenarios=Sc)
  ################################################################################################################### EXPLORE
  
  ################
  ExtractOnly = F
  ################
  
  ### temp
  
  #Sc='D' # one other scenario - the output querying doesn't work with only A
  d=0.3
  #pSVS=0.05 # the proportion of sites receiving a single visit
  
  # Set working directory to where data is to be stored
  dir.create(path ='Output', showWarnings = FALSE)
  
  # Add frescalo weights file to the frescalo directory
  sparta_dir <- paste(find.package('sparta'),'/exec',sep='')
  wts_copy <- file.copy('simWts.txt',sparta_dir)
  if(!wts_copy) stop('Weights file did not copy correctly')
  
  setwd('Output')
  
  xx<-system.time({
    for(d in c(0, 0.3)){ 
      for(pSVS in c(0.05, 0.1, 0.07)) 
      {
        # prepare the run
        ch <- ifelse(d==0, 'V', 'P')
        code = paste(ch,'_',pSVS*100,'SVS_',datecode, sep='') 
        
        if(nSi == 1000){
          set.seed(153) #for continuity, but only for comparing Dutch & British datasets 
        } else {
          Fr <- F
        }
        
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
    }
  })
  
  print(xx)
  
  # Return to homedir
  setwd('..')
  # end
}