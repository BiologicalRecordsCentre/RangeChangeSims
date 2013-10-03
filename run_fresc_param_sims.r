#------------------------------------------------------------------------------
#File:       	run_fresc_param_sims.r
#Description:	R funciton to take simulated occurance data, and return the 
#				time factor and	standard deviation of the time factor for the
#				species 'focal'. Data is read in and moved tot he frescalo
#				folder, frescalo is run and the output file is read in and the 
#				specific details required returned. This code is based on Colin
#				Harrower's run_fresc_param_BRC code. 
#Created:     	06/11/2012
#By:          	Tom August
#Modified:    	
#Version:     	0.1
#------------------------------------------------------------------------------

# REQUIREMENTS
	# Save current directory
		org_wd = getwd()
		
# VARIABLES
	
# FUNCTIONS


# SCRIPT

run_fresc_sims = function(
	fres_in,
	output_dir = "C:/Documents and Settings/TOMAUG/Desktop/Work/BRC/Frescalo - Script/CH_frescalo/Sim",
	frescalo_path = "C:/Documents and Settings/TOMAUG/Desktop/Work/BRC/Frescalo - Script/Frescalo using input file/Frescalo_public_rev2_param_file/Frescalo_2a.exe",
	fres_f_foc = "FocDist.txt", # File name to give the focal data for the frescalo analysis
	fres_f_wts = "simWts.txt", # file name for the frescalo weighting file (needs to already be in directory with frescalo exe file)
	fres_f_log = "Log.txt", # file name for log file created by frescalo (cannot exist in folder already or frescalo will crash)
	fres_f_stat = "Stats.csv", # file name for stats file created by frescalo (cannot exist in folder already or frescalo will crash)
	fres_f_freq = "Freq.csv", # file name for frequencies file created by frescalo (cannot exist in folder already or frescalo will crash)
	fres_f_trend = "Trend.csv", # file name for trend file created by frescalo (cannot exist in folder already or frescalo will crash)
	fres_f_filter = NULL, # file containing sites to be filtered (if required)
	fres_f_nobench = NULL, # file name containing list of species codes of species to exclude from benchmarking (if required needs to already be in directory with frescalo exe file)
	fres_phi_val = NULL, # Parameter value for phi val (NULL equals default of 0.74)
	fres_bench_val = NULL, # Paramter value for local benchmarking threshold (NULL equals default of 0.27)
	spp_names_file = NULL # Path to csv file containing species names in 2 column format, 1) "SPECIES", which contains codes, and 2) "NAME", name belonging to code (if null will assume standard BRC concepts and look them up directly)
	
){
# BODY OF FUNCTION
	# Print status to screen 
		#cat("\nINTIAL SETUP\n",rep("*",20),"\n\n", sep="")
	# Check that frescalo_path exists
		if(!file.exists(frescalo_path)){
			stop("ERROR: supplied frescalo path is not valid")
		}

	# Check output directory exists
		if(!file.exists(output_dir)){
			# If full directory given does not exist but only last folder is missing then create missing directory
			if(file.exists(dirname(output_dir))){
				dir.create(output_dir)
			} else {
				stop("ERROR: output directory path is not valid")
			}
		}
		
	# Create folder to store frescalo data
		fres_dir = file.path(output_dir, "Frescalo")
		if(!file.exists(fres_dir)){
			dir.create(fres_dir)
		}
	
		# Within folder setup input, output and maps_results folders
			# Build paths for directories
			fres_sub_dir = file.path(fres_dir, "Input")
			names(fres_sub_dir) = "INPUT"
			# Create all/missing directories (hiding output from function)
				invisible(sapply(fres_sub_dir[!file.exists(fres_sub_dir)], dir.create))
				
	# Take in data from sims
		filename=file.path(fres_sub_dir["INPUT"], "FocDist.txt")
		write.table(fres_in, file=filename, sep=" ", row.names=FALSE, col.names = FALSE, append = FALSE, quote=FALSE)
	
	# Copy Frescalo file created above to folder where frescalo.exe is stored
		# Extract folder path for the folder in which the frescalo exe is stored
			exe_dir = dirname(frescalo_path)
		# Copy the frescalo file created above into the directory containing the exe
			invisible(file.copy(from = file.path(fres_sub_dir["INPUT"], "FocDist.txt"), to = file.path(exe_dir, "FocDist.txt"), overwrite = TRUE ))
		# Determine if output files from frescalo exist in exe dir 
			f_out_exists = file.exists(file.path(exe_dir, c(fres_f_log, fres_f_stat, fres_f_freq, fres_f_trend)))
			# If so then delete
			if(any(f_out_exists)){
				invisible(file.remove(file.path(exe_dir, c(fres_f_log, fres_f_stat, fres_f_freq, fres_f_trend))[f_out_exists]))
			}
		
	# Create frescalo parameter file (by line)
		# 1) Log fie, 2) data file, 3) weights file, 4) ?no bench file, 5) Site Filter, 6), Stats out file, 7) Freq out file, 8) Trend out file, 9) Phi value, 10) Benchmark value
		# Log file
		cat(fres_f_log,"\n", file=file.path(exe_dir, "params.txt"), append = FALSE)
		# Data File
		cat(fres_f_foc,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		# Weigths
		cat(fres_f_wts,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		# No bench
		cat(fres_f_nobench,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		# Filter
		cat(fres_f_filter,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		# Stats
		cat(fres_f_stat,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		# Freq
		cat(fres_f_freq,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		# Trend
		cat(fres_f_trend,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		# Phi val
		cat(fres_phi_val,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		# Benchmark val
		cat(fres_bench_val,"\n", file=file.path(exe_dir, "params.txt"), append = TRUE)
		
		
	# Run Frescalo using parameter file
		# Print progress
			#cat("\nRUNNING FRESCALO\n",rep("*",20),"\n\n", sep="")
		# Create batch file to change directory and then call frescalo
		if (grepl("linux", R.version$platform)){
			setwd(dirname(frescalo_path))
			system(paste('echo | "',frescalo_path,'"',sep=''))
			setwd(org_wd) 
		    }else{
			cat( paste("cd", normalizePath(dirname(frescalo_path))), "\n", file = file.path(exe_dir, "wincomm.cmd"), append = FALSE)
			cat( basename(frescalo_path), "\n", file = file.path(exe_dir, "wincomm.cmd"), append = TRUE)
			# Run batch file through shell
      		setwd(dirname(frescalo_path))
			shell(shQuote(file.path(exe_dir, "wincomm.cmd")), intern = FALSE)
			setwd(org_wd) 
		    }

    
		# Check log file from frescalo and determine if value from Phi was too low if so then run again with value of Phi + 1 (only if fres_phi_val is NULL)
		if(is.null(fres_phi_val)){
			# Setup inital values of loop variables to ensure entry into while loop (NOTE: only used to enter loop!)
				fres_warn = TRUE
				act_phi = 1.0
				tar_phi = 0.74
			# Setup up while loop to keep going through as long as frescalo gives warning and target is less or equal to than actual value of Phi
			while(tar_phi <= act_phi & fres_warn){
				# Read log file using readLines
				log_out = readLines(file.path(exe_dir, fres_f_log))
				
				# Look for warning in log file
					fres_warn = any(grepl(" [*]{3}[ ]BEWARE[ ][*]{3} ",log_out))
				# Find row stating actual value of phi
					act_txt = log_out[grep("(98.5 percentile of input phi[ ]*)([[:digit:]]{1}[.][[:digit:]]{2})",log_out)]
				# Find row stating target value of phi
					tar_txt = log_out[grep("(Target value of phi[ ]*)([[:digit:]]{1}[.][[:digit:]]{2})", log_out)]
				# Read value of Phi
					act_phi = as.numeric(gsub("(98.5 percentile of input phi[ ]*)([[:digit:]]{1}[.][[:digit:]]{2})","\\2", act_txt))
				# Read Target value of phi
					tar_phi = as.numeric(gsub("(Target value of phi[ ]*)([[:digit:]]{1}[.][[:digit:]]{2})","\\2", tar_txt))
					
				if(tar_phi <= act_phi & fres_warn){
					# Print notice to screen to say will need to rerun
						cat("NOTE: Targest value of phi may be too small, frescalo will be rerun with a larger target value of Phi\n\n")
					# Read in params file and alter 9 value (Phi value line)
						param_temp = readLines(file.path(exe_dir,"params.txt"))
					# Set target phi to act_phi + 0.01 to make sure rounding issues don't cause failure again
						param_temp[9] = act_phi + 0.01
					# Write lines to param file
						writeLines(param_temp, file.path(exe_dir,"params.txt"))
						
					# Determine if output files from frescalo exist in exe dir 
						f_out_exists = file.exists(file.path(exe_dir, c(fres_f_log, fres_f_stat, fres_f_freq, fres_f_trend)))
						# If so then delete
						if(any(f_out_exists)){
							invisible(file.remove(file.path(exe_dir, c(fres_f_log, fres_f_stat, fres_f_freq, fres_f_trend))[f_out_exists]))
						}
						
					# Rerun frescalo
						# Create batch file to change directory and then call frescalo
						if (grepl("linux", R.version$platform)){
						setwd(dirname(frescalo_path))
						system(frescalo_path)
						setwd(org_wd) 
						}else{
						cat( paste("cd", normalizePath(dirname(frescalo_path))), "\n", file = file.path(exe_dir, "wincomm.cmd"), append = FALSE)
						cat( basename(frescalo_path), "\n", file = file.path(exe_dir, "wincomm.cmd"), append = TRUE)
						# Run batch file through shell
						setwd(dirname(frescalo_path))
						shell(shQuote(file.path(exe_dir, "wincomm.cmd")), intern = FALSE)
						setwd(org_wd) 
						}
				}
			}
		}
			
	
	
	# Extract data for focal species
	
			# Set file path for trends output
			trend_path = file.path(exe_dir, fres_f_trend)
			
			# Read file
		
		    # Read majority of data (skip 1st line)
		    trend <<- read.csv(trend_path, header=TRUE, stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA", "*******","********","1.#QO","-1.$","-1.#IO","1.#R","1#INF"), comment.char = "")
	    	
			# Extract time factors and standard errors for the focal species
			focal_rows<-subset(trend,Species__=='focal',select=c(Time______,TFactor,St_Dev))
			# Create a dataframe of time period, Time factors and St_Dev. Sorted by time period
			Tfac_StDev<-focal_rows[with(focal_rows, order(Time______)), ]
		    
        colnames(Tfac_StDev)<-c('TimePeriod','Tfactor','StDev')
	
	# NI: return the actual value of Phi
	attr(Tfac_StDev, 'Phi') <- act_phi
	
	# Return vector of Time factor and St_Dev of focal species
		return(Tfac_StDev)
}
