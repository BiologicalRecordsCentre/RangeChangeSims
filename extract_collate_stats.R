
# Nick Isaac 24 January - 5 February 2013
# Extract and explore the simulation results
# the simulations are now very slow
# so I'm running each set of input values for only a small number of reps and building up gradually
# these are the functions to extract results from multiple sets

# 3 May: rewrote the extract_output function (didn't work previously)

# NEED TO TO:
# Add number of reps somewhere (table?)
# make it prettier
# add a separate version for querying output, not stats. Then stop writing csv stats files.
# the output has key parameters included as attributes: it would be more robust to use these.

# This version of the function is for Stats files produced before and after 4/2/13
# On that date I changed the naming structure, swapping the position of the date code (e.g. 130204) with the test id ('V' or 'P')
# This works on both versions because the positions of nVisits, nSpp & nReps are unchanged


extract_output <- function(pathname='.'){
    # read the output files directly
    files <- list.files(path=pathname, pattern='SimOutput') # all the rData files with error rates
    output <- lapply(files, function(filename) {
        load(paste(pathname,filename,sep='/')) # -> output
        
        # convert the output into a single long dataframe
        moutput <- melt(output)
        names(moutput)[1:3] <- c('Statistic','Scenario','repnum') #4th column remains called 'value'
        
        # now put all the informatin about the simulation parameters into s
        moutput <- (cbind(moutput, t(unlist(attr(output, "simpars")))))
        moutput$id <- attr(output, "id")
        
        return(moutput)
    })
    cols <- 1:ncol(output[[1]]) # each element of the list should have the same number of columns
    return(melt(output, id=cols)[,cols])
}

extract_error_rates <- function(pathname='.', chr='V'){
    require(reshape2)
    # find all the files with error rates
    files <- list.files(path=pathname, pattern='SimStats') # all the csv files with error rates
    files <- files[grepl(paste('_',chr,'_',sep=''), files)] # restrict those to Type I error rates (not power test)
    
    rates <- sapply(files, USE.NAMES=F, simplify=F, function(filename) {
        # read in the file of error rates and format the error rates in a sensible format
        stats <- read.csv(paste(pathname,filename,sep='/'), stringsAsFactors = FALSE)
        error_rates <- stats[grep('plr', stats$X) : grep('^nSites$', stats$X),] # extract just the error rates (removing the summary stats)
        error_rates <- melt(error_rates) # convert to long format
        
        # now extract simulation parameters from the filename
        underscores <- gregexpr('_', filename)[[1]]
        error_rates$nVisits <- as.numeric(substr(filename, underscores[2]+1, underscores[3]-2))
        error_rates$nSpp <- as.numeric(substr(filename, underscores[3]+1, underscores[4]-3))
        error_rates$nReps <- as.numeric(substr(filename, underscores[5]+1, regexpr('\\.', filename)[[1]]-1))
        return(error_rates)
    })
    x <- melt(rates, id=c(1,2,4,5,6))
    names(x)[1:2] <- c('Method','Scenario')
    x$numsig <- with(x, nReps * value)
    
    # just in case there are multiple runs with the same parameter values, we'll aggregate
    a <- aggregate(numsig ~ Method + Scenario + nVisits + nSpp, data=x, sum)
    b <- aggregate(nReps ~ Method + Scenario + nVisits + nSpp, data=x, sum)
    rates <- merge(a, b)
    rates$error_rate <- with(rates, numsig/nReps)
    return(rates)
}


extract_recording_rates <- function(pathname='.', chr='V'){
    require(reshape2)
    # find all the files with error rates
    files <- list.files(path=pathname, pattern='SimStats') # all the csv files with error rates
    files <- files[grepl(paste('_',chr,'_',sep=''), files)] # restrict those to Type I error rates (not power test)
    
    rates <- sapply(files, USE.NAMES=F, simplify=F, function(filename) {
        # read in the file of error rates and format the error rates in a sensible format
        stats <- read.csv(paste(pathname,filename,sep='/'), stringsAsFactors = FALSE)
        rates <- stats[grep('Recs_total', stats$X),] # extract just the error rates (removing the summary stats)
        #rates <- melt(rates) # convert to long format
        
        # now extract simulation parameters from the filename
        underscores <- gregexpr('_', filename)[[1]]
        rates$nVisits <- as.numeric(substr(filename, underscores[2]+1, underscores[3]-2))
        rates$nSpp <- as.numeric(substr(filename, underscores[3]+1, underscores[4]-3))
        rates$nReps <- as.numeric(substr(filename, underscores[5]+1, regexpr('\\.', filename)[[1]]-1))
        return(rates[,-1])
    })
    x <- melt(rates, id=c('nVisits','nSpp','nReps'))
    x$sumrecs <- with(x, value * nReps)
    names(x)[4] <- c('Scenario')
    
    # just in case there are multiple runs with the same parameter values, we'll aggregate
    a <- aggregate(sumrecs ~ Scenario + nVisits + nSpp, data=x, sum)
    b <- aggregate(nReps ~ Scenario + nVisits + nSpp, data=x, sum)
    rates <- merge(a, b)
    rates$Recs_rep <- with(rates, sumrecs/nReps)
    rates$Recs_spp_yr <- with(rates, Recs_rep/(nSpp*10))
    
    return(rates)
}
