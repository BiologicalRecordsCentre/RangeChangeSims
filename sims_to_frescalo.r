# Tom August - 06/11/12

# Function to take output from simualtion functions (a dataframe)
# and turn them into a format suitable for Frescalo i.e.
# three columns: Location | Species | Time Period

# NOTE: The weights file must be placed in the directory with Frescalo_2a.exe
# and named Wts.txt

# records is an element of a list (also named records) produced by
# iterate_all_scenarios() in 'Sim functions'
# output_dir is the directory where a text file is temporarily stored
# before running Frescalo
# frescalo_path is the full path to Frescalo_2a.exe
# tp2, if true, converts bins the time periods into two periods about the 
# mid year in the series


sims_to_frescalo <- function (records,output_dir,frescalo_path,tp2=FALSE){

  # take columns Frescalo needs
  cols<-c("Site","Species","Year")
  fres_in<-records[cols]
  
  # bin years into two time periods
  if(tp2){
      split=mean(unique(fres_in$Year))
      lower<-fres_in$Year<split
      upper<-fres_in$Year>split
      fres_in[lower,]$Year<-1
      fres_in[upper,]$Year<-2
      }
  
  # call the frecalo function
  Tfac_StDev<-run_fresc_sims(fres_in,output_dir,frescalo_path)
  
  # return the time factor and standard deviation for the species 'focal'
  return(Tfac_StDev)
  }