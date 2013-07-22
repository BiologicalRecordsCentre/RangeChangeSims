# Nick Isaac
# 22 Ocobter 2012
# 
# Complete rewrite of most simulation functionality
# 
#################################################
#
# TO DO
# Check plot functions with multiple MM thresholds
#
#################################################
require(reshape2)
require(lme4)

#Source Frescalo code (must be in a subdirectory of the wd)
source('Frescalo/Process files/run_fresc_param_sims.r')
#source('Frescalo/Process files/sims_to_frescalo.r') # moved further down

#################################################
#################################################
################################################# GENERIC FUNCTIONS - TO BE ADDED TO TO 'RANGE CHANGE FUNCS'
occurrence <- function(x) length(x) > 0 # takes a vector and returns whether the length is greater than 0


cast_recs <- function(records, resolution='visit', focalspname='focal'){
	#takes a set of records and returns a dataframe suitable for analysis
	#columns include the list length and presence/absence of the focal species
    #6 December: minor alterations to make it work with nonfocal species

	if(resolution=='visit'){
		castrecs <- dcast(records, Year + Site + Visit ~ ., value.var='Species', fun=LenUniq)
		spdata <- subset(records, Species==focalspname) #simply the subset of data for the focal species

    } else if(resolution=='kmyr'){
		castrecs <- dcast(records, Year + Site ~ ., value.var='Species', fun=LenUniq)
		spdata <- subset(records, Species==focalspname) #simply the subset of data for the focal species
		spdata <- dcast(spdata, Year + Site ~ ., value.var='Species', fun=occurrence)
	}
	names(castrecs)[ncol(castrecs)] <- 'L'
	# merge the two
	castrecs <- merge(castrecs, spdata, all=T, nomatch=FALSE)
 
    names(castrecs)[ncol(castrecs)] <- focalspname #changed from 'focal'
	castrecs[ncol(castrecs)] <- as.logical(!is.na(castrecs[ncol(castrecs)]))
    return(castrecs)
}


frescalo_bootstrap <- function(fr, n=1000){
  #takes the output from sims_to_frescalo()
  #it returns a trend and p-value, estimated from bootstrapping
  #NB 'trend' estimates are measured in units of 'time periods', not the underlying years (ok if nTp=nYr)
  require(reshape2)
  x <- apply(fr, 1, function(x) rnorm(n, mean=x[2], sd=x[3]))
  mod <- summary(lm(value ~ Var2, data=melt(x)))
  return(c(mod$coefficients['Var2','Estimate'], mod$coefficients['Var2','Pr(>|t|)']))
}


frescalo_trend<- function(fr){
  #takes the output from sims_to_frescalo()
  #it returns a trend and p-value
  #if nTP=2 then it uses a z-test
  #if nTP>2 is calculates a linear trend, ignoring stdev
  #NB 'trend' estimates are measured in units of 'time periods', not the underlying years (ok if nTp=nYr)
  if(nrow(fr)==2){
    trend <- diff(fr$Tfactor)
    z <- trend/sqrt(sum(fr$StDev^2))
    p <- one_to_two_tail(pnorm(z))
  } else {
    mod <- summary(lm(Tfactor ~ TimePeriod, data=fr))
    trend <- mod$coefficients['TimePeriod','Estimate']
    p <- mod$coefficients['TimePeriod','Pr(>|t|)']
  }
  return(c(trend, p))
}


Convert_records_to_2tp <- function(records, splityr) {
	# takes raw biological records and generates a dataframe suitable for Telfer etc
	n1 <- with(subset(records, Year <= splityr), tapply(Site, Species, LenUniq)) # number of sites for each species in t1
	n2 <- with(subset(records, Year > splityr), tapply(Site, Species, LenUniq))
	d1 <- LenUniq(subset(records, Year <= splityr)$Site) # total number of sites in t1
	d2 <- LenUniq(subset(records, Year > splityr)$Site) # total number of sites in t2
	gridcell_counts <- data.frame(n1,n2)
	gridcell_counts[is.na(gridcell_counts)] <- 0 # fix the bug of unrecorded species (added 10/1/13)
    attr(gridcell_counts, 'denom') <- c(d1,d2)
	return(gridcell_counts)
}


fit_LadybirdMM <- function(MMdata, nsp=2, nyr=3){ #0.44 seconds
	# this version returns only the coefficients and the proportion of observations that were 'well-sampled'
	# 14 December: minor change: I removed the call to MMdata to outside the function.
    #   this Gives extra flexibility to analyse disaggregated data.
    require(lme4)

    #subset the data
    i <- MMdata$L >= nsp
	i[i==T] <- is.gridcell.wellsampled(MMdata$Site[i], n=nyr)

    # 4/2/13: with small numbers of visits, it's possible that no gridcells meet the criterion
        # I also centred the data on the mid year in order to improve numerical convergence 
    my <-median(unique(MMdata$Year))
    if(sum(i) > nyr){ # PERHAPS CHANGE THIS TO LenUniq(MMdata[i,]$Site) > 6
        x <- try({# another source of the bug error is if there's not enough to fit the model
            MM <- glmer(focal ~ I(Year-my) + (1|Site), MMdata, subset=i, family=binomial)
            coefs <- as.numeric(summary(MM)@coefs[2,])
            }, silent=T)
        if(class(x)=='try-error') save(list(MMdata, nsp), file='MM_ber_tryerror.rData')
        # end bug check
    } else coefs <- rep(NA, 4)
    
	return(c(coefs, sum(i)/length(i))) # keep this a two step process in case we later decide to extract other info
}


fit_ladybirdMM_bin <- function(indata, nsp=2, nyr=3, od=F, V=F){
    # indata (simdata) contains 1 row per valid visit, with a true or false for whether it was recorded
    # now we aggregate this to the year:monad level, with a binomial response (not bernoulli) <- thanks to Ben Bolker
    # initially I calculated nVR directly from records, but this includes all the lists of length 1 and the poorly-sampled grid cells
    # optional arguments for modelling Overdispersion (very slow) and verbose output

    require(lme4)
    
    #subset the data: remove the short lists (defined by nsp)
    data <- subset(indata, L>=nsp)
    
    # of these visits, which are on well-sampled sites?
    # 21 June: if no thresholds are applied then no need to look for well-sampled set
    if(nsp>1 & nyr>1) data <- subset(data, is.gridcell.wellsampled2(data, n=nyr))
    
    if(nrow(data) > 0){
        MMdata <- dcast(data, Year + Site ~ ., fun=length, value.var='L') #how many lists?
        names(MMdata)[ncol(MMdata)] <- 'nVis'
    
        MMdata$nVR <- as.numeric(acast(data, Year + Site ~ ., fun=sum, value.var='focal'))
    
        # to fit a binomial model, we define a new column containing the number of visits without an observaiton of the focal
        MMdata$failures <- with(MMdata, nVis - nVR)
    
        # centre the Year on the median value (for numerical stability)
        MMdata$cYr <- MMdata$Year - median(unique(MMdata$Year))

        x <- try({# another source of the bug error is if there's not enough to fit the model
            if(od) {
                MMdata$obs <- 1:nrow(MMdata)
                MM <- glmer(cbind(nVR, failures) ~ cYr + (1|Site) + (1|obs), data=MMdata, family=binomial, verbose=V)
            } else {
                MM <- glmer(cbind(nVR, failures) ~ cYr + (1|Site), data=MMdata, family=binomial, verbose=V)
            }    
            coefs <- as.numeric(summary(MM)@coefs[2,])
            }, silent=T)
        if(class(x)=='try-error') save(MMdata, file='MM_bin_tryerror.rData')
        # end bug check

        # calculate the number of rows in the model as a proportion of the total site-year combos
        p_used <- nrow(MMdata)/nrow(unique(indata[,c('Year','Site')]))
        result <- c(coefs, p_used)
    } else result <- c(rep(NA, 4),0) # there's no 'well-sampled' cells
    
    return(result) # keep this a two step process in case we later decide to extract other info
}

fit_Maes <-function(records, splityr, min_sp=5){
    # fits the method described in Maes et al 2012 (Biol Cons 145: 258-266)
    # I checked this by comparing the results against those in Maes et al SI
    # I've added a sampling theory to estimate p-values. We'll see whether it's robust!
    # the function returns numbers of sites in 2 time periods (out of the well-sampled ones), the %trend and p-value
    require(reshape2)
    
    # which time period is each record in?
    records$tp <- 1+ (records$Year > splityr)
    
    # convert the records into a 3D array
    rc <- acast(records, Species ~ Site ~ tp, fun=occurrence, value.var=2)
    
    # what is the observed species richness in each cell in each year
    rc1 <- apply(rc, c(2,3), sum)
    
    # the number of sites in each time period
    nS <- colSums(rc1>0)
    
    # which sites have are well-sampled? (defined as having at least min_sp species in BOTH time periods)
    well_sampled <- as.numeric(dimnames(rc1)[[1]][apply(rc1, 1, function(x) all(x>=min_sp))])
    
    # look at just the data for these well-sampled cells
    rc2 <- rc[,dimnames(rc)[[2]] %in% well_sampled,]
    
    # how many sites for each species in each time period?
    rc3 <- apply(rc2, c(1,3), sum)
    
    # calculate the relative distribution in each time period
    rd1 <- rc3[,1]/nS[1]
    rd2 <- rc3[,2]/nS[2]
    
    trend <- 100 * (rd2-rd1)/rd1
    
    # we can assess the significance of this method as follows
    # first we assume the distribution of each species is equal among poorly-sampled and well-sampled sites
    # Under the null hypothesis that total proportion of sites has not changed,
    # we can calculate the binomial probability of the 'estimated number of successes
    # estimated number of successes is defined here as nsr2 (number of sites recorded in t2)
    # where the number of trials = nS[2]
    
    nsr2 <- nS[2] * rc3[,2] / length(well_sampled)
    true_probs <- rc3[,1]/length(well_sampled)
    
    pval <- mapply(FUN=pbinom, q=nsr2, prob=true_probs, MoreArgs=list(size=nS[2]))
    
    #these are one-tailed: convert them to one-tailed
    pval <- one_to_two_tail(pval)
    
    Maes <- data.frame(N1=rc3[,1], N2=rc3[,2], trend=trend, pval=pval) 
    attr(Maes, 'nSites') <- nS
    attr(Maes, 'wellsampled') <- length(well_sampled)
    return(Maes)
}

#takes a dataframe of observations, with a vector indexing the time period
fit_Telfer <- function(gridcell_counts, min_sq=5) {
	#takes output from Convert_to_2tp and fits the telfer model
	#convert to proportions	
	p1 <- gridcell_counts$n1/attr(gridcell_counts, 'denom')[1]
	p2 <- gridcell_counts$n2/attr(gridcell_counts, 'denom')[2]

	#Telfer's method is the *standardized* residual from a logit-logit regression
	model1 <- lm( log(p2/(1-p2)) ~ log(p1/(1-p1)), subset=gridcell_counts$n1 >=min_sq & p2 >0)
	return(rstandard(model1))
	}

is.gridcell.wellsampled <- function(CellID, n=3){
	# modified version of filter.gridcells()
	#takes the dataframe and returns the rownumbers of the dataset identifying well-sampled the gridcells
	nyrs <- table(CellID)
	return(CellID %in% names(nyrs)[nyrs >= n])
	}

is.gridcell.wellsampled2 <- function(data, n=3){
    #takes the dataframe and returns a logical vector identifying well-sampled gridcells
    #this version copied from the Odonata trend analysis (late Feb 2013)
    #changed 4 March to filter on number of years, not number of visits
    require(reshape2)
    x <- acast(data, Site~Year, fun=length, value.var='L') # x is the number of visits
    num_yrs_visits  <- rowSums(x>0) # convert nVisits to binary and sum across years 
    sites_to_include <- num_yrs_visits >= n  
    return(data$Site %in% dimnames(x)[[1]][sites_to_include])
}
################################################# SIM-SPECIFIC FUNCTIONS

bias_focal <- function(records, disadvantage=0.5){
	# a clever way to simulate bias *against* is to delete records from 'well-sampled' set by stratified random sampling
	# disadvantage defines the proportion of lists where the focal species should be 'unrecorded'
	# If a single number then it's constant across years.
	# if numeric, this remains constant over time. If two numbers, the second number identifies the FINAL value
	increment <- ifelse(length(disadvantage)==1, 0, diff(disadvantage) / max(records$Year)) # which extra ones should be included in each year
	records$record_id <- 1:nrow(records)
	focal <- subset(records, subset=Species=='focal')
	p_lost <- disadvantage[1] + increment * focal$Year # probability of being 'lost' (un-recorded) is a function of year and disadvantage
	lost <- rbinom(n=nrow(focal), size=1, prob=p_lost) # a vector of zero & 1s indicating  
	record_id_to_lose <- focal$record_id[lost==1]
	return(records[-record_id_to_lose,]) # 
	}


create_data_1 <- function(nSites=500, nSpecies=25, pFocal=0.5){
	# The original version
    # 24 October: I removed the 'ubiquitous' option & replaced it with a normal distribution of site occupancy
	# creates a site by species matrix of presence-absence data
	#focal <- sample(x=c(1,0), size=nSites, replace=T, prob=c(pFocal, 1-pFocal))
	focal <- rbinom(n=nSites, size=1, prob=pFocal)
	#if(ubiquitous) x <- matrix(data=1, nr=nSites, ncol=nSpecies) #all non-focal species are everywhere
	#else x <- matrix(data=sample(x=c(1,0), size=nSites*nSpecies, replace=T), nr=nSites, ncol=nSpecies) #all nonfocal species have 50% probability of occupancy
	#else x <- matrix(data=rbinom(n=nSites*nSpecies, size=1, prob=0.5), nr=nSites, ncol=nSpecies) #all nonfocal species have 50% probability of occupancy
	
	# NEW: the probability of occupancy for any species is changed from a constant (0.5) to a random normal deviate
	prob <- rnorm(n=nSpecies, mean=0.5, sd=0.25)
	prob[prob>=1] <- 1
	prob[prob<0.05] <- 0.05

	x <- t(replicate(nSites, rbinom(n=nSpecies, size=1, prob=prob)))
	#summary(colSums(x)/nSites) # matches well with prob
	dimnames(x)[[1]] <- paste('site',1:nSites, sep='')
	dimnames(x)[[2]] <- paste('spp',1:nSpecies, sep='')
	attr(x,'richness') <- summary(rowSums(x))
	x <- cbind(focal, x)
}


create_data_2 <- function(nSites=500, nSpecies=25, pFocal=0.5, bsh1=6, bsh2=4, k=3){
  # 3 December: as 'old' above, except
  # a) nonfocal species are defined by a beta distribution (occupancy & detectability)
  # b) includes an attribute of the true detection probability for each species
    
  focal <- rbinom(n=nSites, size=1, prob=pFocal)
    
  # define the probability of occupancy based on a beta distribution
  p_detect <- rbeta(n=nSpecies-1, bsh1, bsh2)
  #sum(p_detect[p_detect>=quantile(p_detect, 1-0.27)])/sum(p_detect) # around 0.34, ie 34% of detectablity lies within benchmarks. target is 45%
  p_occ <- p_detect ^ k
    
  x <- t(replicate(nSites, rbinom(n=nSpecies, size=1, prob=p_occ)))
  dimnames(x)[[1]] <- paste('site',1:nSites, sep='')
  dimnames(x)[[2]] <- paste('spp',1:nSpecies, sep='')
  x <- cbind(focal, x)
  attr(x,'richness') <- as.numeric(rowSums(x))
  attr(x,'p_detect') <- c(pFocal, p_detect)
  return(x)  
}


create_data_3 <- function(nSites=500, nSpecies=25, pFocal=0.5, bsh1=6, bsh2=4, k=3){
    # 3 December: as 'old' above, except
    # a) nonfocal species are defined by a beta distribution (occupancy & detectability)
    # b) includes an attribute of the true detection probability for each species
    # 10 January: Further modified so that detectablity follows Mark Hill's formula for detectability
    # I also removed the correlation between abundance and occupancy 
    
    focal <- rbinom(n=nSites, size=1, prob=pFocal)
    
    # define the probability of occupancy based on Mark Hill's formula
    Hill_detectability <- function(nSp, a=2.005, b=-2.545){
        # Here I apply a global version of Mark Hill's empirically-derived detectability function
        # It's appropriate for neighbourhoods where the expected number of species is 1/2 the total species pool
        # f’(R’) = 1 - exp(-exp(a+b R’))
        # Parameters a and b were derived empirically from Bryophyte data
        # R is simply species number divided by 1/2 the total species richness
        # this way, the median detectability is about 0.42 and the minimum around 0.045 (using Mark's a and b)
        R <- 1:nSp/(0.5*nSp)    
        pd <- 1 - exp(-exp(a+b*R))
        #sum(pd[pd>=quantile(pd, 1-0.27)])/sum(pd) # around 0.34, ie 34% of detectablity lies within benchmarks. target is 45%
    }
    
    p_detect <- Hill_detectability(nSpecies-1) # substract focal since we're fixing that constant
    #p_occ <- p_detect ^ k
    p_occ <- rbeta(n=nSpecies-1, bsh1, bsh2) ^ k
    
    x <- t(replicate(nSites, rbinom(n=nSpecies, size=1, prob=p_occ)))
    dimnames(x)[[1]] <- paste('site',1:nSites, sep='')
    dimnames(x)[[2]] <- paste('spp',1:nSpecies, sep='')
    x <- cbind(focal, x)
    attr(x,'richness') <- as.numeric(rowSums(x))
    attr(x,'p_detect') <- c(pFocal, p_detect)
    return(x)  
}


create_data <- function(nSites=500, nSpecies=25, pFocal=list(Occ=0.5, DetP=0.5), ONR=F, bsh1=2, bsh2=2){
    # 3 December: as 'old' above, except
    # a) nonfocal species are defined by a beta distribution (occupancy & detectability)
    # b) includes an attribute of the true detection probability for each species
    # 24 January: I fixed some inconsistency in the use of nSpecies
    #   For simplicity, the User-defined number of species is the number of NON-FOCALS
    #   I've also changed the shape of the interspecific occupancy distribution to be symmetrical around 0.5
    #   Although [mean occupancy]<<0.5 for most groups, this high number is necessary to generate enough 'global benchmarks'
    #   (i.e. species that are common enough everywhere to give Frescalo a decent change of generating neighbourhoods)
    #   And added an optional ONR (Occupancy Abundance(Number) Relationship) - the N follows McGill
    #   When False (the default) there is NO correlation between occupancy and detectability (abundance)
    #   When True, there is a perfect correlation between Detectability and Occupancy
    # 22 January: I further modified the detectability function, following consultation with Mark Hill.
    #   The function descrbed in create_data_3 (the last current version) is for a survey of intensity=1
    #   For bryophytes, surveying a hectad to this intensity takes 8 days!
    #   In Mark's concept, the pool of possible species being detected at a site is infinite,
    #       but the 'expected number of species from a complete survey (intensity >> 1) is known
    #   I want to simulate a standard visit, which is much less intense, but where the most detectable species are highly likely to be recorded
    #   I tried using an exponential decay function but it was much too steep.  
    # 9 April: pFocal separated into Occupancy & probability of detection
    
    # first define the vector of presence-absence for the focal species
    focal <- rbinom(n=nSites, size=1, prob=pFocal$Occ)
    
    HillP <- function(nSp, pse=1, a=2.005, b=-2.545, S=0.3){
        # generates a detectability decay curve as per Mark hill's Frescalo paper
        # under the default settings, the sppSVS has an 89% detection prob, the least detectable 16%.
        # the median detectability is 46%
        # this decay curve is close to linear
        R <- 1:nSp/(pse*nSp)  
        pd <- 1-exp(-S*exp(a+b*R))
    }
    
    p_detect <- HillP(nSpecies)
    if(ONR) p_occ <- p_detect
    else p_occ <- rbeta(n=nSpecies, bsh1, bsh2)
    
    x <- t(replicate(nSites, rbinom(n=nSpecies, size=1, prob=p_occ)))
    dimnames(x)[[1]] <- paste('site',1:nSites, sep='')
    dimnames(x)[[2]] <- paste('spp',1:nSpecies, sep='')
    x <- cbind(focal, x)
    attr(x,'richness') <- as.numeric(rowSums(x))
    attr(x,'p_detect') <- c(as.numeric(pFocal$DetP), p_detect)
    return(x)  
}


generate_all_scenarios <- function(nSites=1000, nSpecies=50, nYrs=10, pSVS=0.05, mv=20, vrs=F, stoch=T, Scenarios='BCDF', combos=F,
                                   save_data=F, id='', pFocal=list(Occ=0.5, DetP=0.5),p_short=list(init=0,final=0.2), pDetMod=0.2, decline=0){
    # This is a wrapper for creating data
    # 25 October: added 'decline' for testing power
    #  3 December: removed B1 & replaced E
    # 25 January 2013: I added B3: new sites are biased toward the focal species
    #  4 February: added F: decline in a nonfocal species
    #  7 February: added the option of user-controlled subset of scenarios
    # 13 February: I added the optional argument pFocal to allow exploration of changing it from outside
    # 14 February: I dropped D1 and F1
    # 8 April: removed B3n
    
    # the actual data   
    true_data <- create_data(nSites=nSites, nSpecies=nSpecies, pFocal=pFocal)
    
    #cntrlsp <- median_occ(true_data) #the identity of a control species
    records_A = generate_records(nYrs=nYrs, pSVS=pSVS, true_data=true_data, decline=decline, mv=mv, vrs=vrs, stoch=stoch)
    records <- list(A_EvenRcrdng=records_A)
    
    if(grepl('B', Scenarios)){ # run the B Scenarios
        #records$B1_IncrnSite = subset_recs_by_X(records_A, X=records_A$Site, init=0.5), # first year has 50% sites, final year has all
        records$B2_IncrnVisit = subset_recs_by_X(records_A, X=records_A$Visit, init=0.5) # first year has 50% visits, final year has all
        records$B3f_IncrnVBiasFc = subset_recs_Biased(records_A, init=0.5, target=true_data[,1])
        #records$B3n_IncrnVBiasNf = subset_recs_Biased(records_A, init=0.5, target=true_data[,cntrlsp]) # bias wrt a species of median detectability      
    }
    if(grepl('C', Scenarios)){ # run the C Scenarios
        records$C1_pShortLEven = shorten_lists(records_A, p_short=mean(unlist(p_short))) #a fixed proportion will be short, split evenly between incidental and others 
        records$C2_pShortLIncr = shorten_lists(records_A, p_short=p_short) #proportion of short lists increases over time
    }
    if(grepl('D', Scenarios)){ # run the D Scenarios
        #records$D1_SelectvEven = bias_focal(records_A, disadvantage=pDetMod) #Remove records of the focal species at random (bias *against* detection) 
        records$D2_SelectvIncr = bias_focal(records_A, disadvantage=c(pDetMod,0)) #Bias against detection decreases to zero (i.e. increasing apparency)
    }  
    if(grepl('E', Scenarios)){ # generally not used - after 25/4 it is unuseable since we've dropped nVisits as a parameter
        #records$E1n_IncidentalNfEv = rbind(records_A, incidental_records(true_data, nYrs=nYrs, nVisits=nVisits, p_short=mean(unlist(p_short)), focal=F)) #add nonfocal records
        #records$E2n_IncidentalNfIn = rbind(records_A, incidental_records(true_data, nYrs=nYrs, nVisits=nVisits, p_short=p_short, focal=F)) #add nonfocal records
        #records$E1f_IncidentalFcEv = rbind(records_A, incidental_records(true_data, nYrs=nYrs, nVisits=nVisits, p_short=mean(unlist(p_short)), focal=T)) #add nonfocal records
        #records$E2f_IncidentalFcIn = rbind(records_A, incidental_records(true_data, nYrs=nYrs, nVisits=nVisits, p_short=p_short, focal=T)) #add nonfocal records
    }
    if(grepl('F', Scenarios)){
        # Added 4/2/13: A new scenario in which a nonfocal species (indexed by cntrlsp) declines dramatically (50% over the period)
        #records$F_NfDecline <- generate_records(nYrs, pSVS, true_data, decline=c(decline,0.5), which.decline=c(1,cntrlsp))
        #records$F_NfDecline10 <- generate_records(nYrs, pSVS, true_data, decline=c(decline,0.5), which.decline=0.1)
        records$F_NfDecline <- generate_records(nYrs=nYrs, pSVS=pSVS, true_data=true_data, decline=c(decline,0.3), which.decline=0.5, mv=mv, vrs=vrs, stoch=stoch) # half the species are VU
    }
    if(combos){# now add the combination scenarios
        records$B2C1_IncVsShtEv = shorten_lists(records$B2_IncrnVisit, p_short=mean(unlist(p_short)))
        records$B2C2_IncVsShtIn = shorten_lists(records$B2_IncrnVisit, p_short=p_short)
        #records$B2D1_IncVsSelEv = bias_focal(records$B2_IncrnVisit, disadvantage=pDetMod)
        #records$C1D1_ShtEvSelEv = bias_focal(records$C1_pShortLEven, disadvantage=pDetMod)
        #records$C2D1_ShtInSelEv = bias_focal(records$C2_pShortLIncr, disadvantage=pDetMod)
        #records$B2C1D1_IncVShtSelEv = bias_focal(records$B2C1_IncVsShtEv,disadvantage=pDetMod) #should kill LL model (new 6/11)
    }
    # for validity, we there's an alternative version of B2, done differently (for reality check). Doesn't work for power test
    #if(decline==0) records$B2_IncrnVisit_a = generate_records(nYrs, pSVS=c(pSVS/2,pSVS/(2*nYrs)), true_data)
    
    attr(records,"true_data") <- true_data
    attr(records,"simpars") <- list(nSpecies=nSpecies,nSites=nSites,nYrs=nYrs,pSVS=pSVS,pFocal=pFocal,p_short=p_short,
                                    pDetMod=pDetMod,decline=decline,visit_rich_sites=vrs, stochastic=stoch, max_vis=mv)
    if(save_data) save_records(records, id)
    
    return(records)
}


generate_records <- function(nYrs=2, true_data, decline=0, which.decline=1, sites_to_use=c(1,0), 
                             pSVS=0.05, mv=20, vrs=F, stoch=T) {
    
    # wrapper for recording_cycle(), allowing it to be run over many years
    # sites_to_use was intended to be the way I implemented scenario B1, but no longer necessary
    # 25 October: added 'decline' for testing power
    # nVisits is allowed to vary among years. If one value is supplied, the visit rate remains constant
    # If two values are supplied then the first value is the initial rate, and the second value is the annual increment
    # 4/2/13: I modified this to allow a species other than the focal to be in decline (which.decline). Default=1 (focal)
    # 5/2/13: Further modification to allow two species to decline simultaneously
    # 6/2/13: Additional statement to allow a variable number of declining species
    #       if which.decline < 1 then it specifies a proportion of species to select at random
    #       otherwise, which.decline specifies the identities of the species to decline
    #       Under scenario F, decline must be supplied as a vector of length 2
    #       The first element specifies the decline rate of the focal species, the second selected nonfocals
    # 24/4/13: modified so that nVisits is replaced by pSVS, the proportion of sites receiving a single visit
    
    if(length(pSVS)==1) pSVS <- c(pSVS[1], 0)
    
    if (sum(decline)==0) {records <- lapply(1:nYrs, function(i) 
        recording_cycle(pSVS=pSVS[1]+i*pSVS[2], true_data=true_data, max_vis=mv, VisRichSites=vrs, stochastic=stoch)) #can start loop at zero
    } else {
        # This part of the loop is triggered under two circumstances
        # 1: Testing power
        # 2: when 1 or more nonfocal (nuisance) species are declining (scenario F)
        # To test power, we simulate a declining population
        # For the IUCN criterion A, a decline of 30% over 10 years leads to a Vulnerable listing.
        
        # When nonfocal species are declining, which.decline specifies which ones
        if(length(which.decline) == 1) if(which.decline < 1) { # double if avoids a warning message
            #which.decline is a proportion of species, to be selected at random
            #convert the proportion into a number
            nsp_declining <- round(which.decline * (ncol(true_data)-1))
            # select that number at random from the nonfocal species
            which.decline <- sample(2:ncol(true_data), size=nsp_declining, replace=F)
            
            # is the focal also declining? (i.e. are we testing power under scenario F?)
            if(decline[1] == 0) # NO
                decline <- rep(decline[2], nsp_declining)
            else {#Yes
                decline <- c(decline[1], rep(decline[2], nsp_declining))
                which.decline <- c(1, which.decline)
            }
        }        
        if(length(which.decline) != length(decline)) stop('incompatible lengths')
        
        # each year, every record has a finite probability of persisting to the next year
        ann_persit_rate <- (1-decline)^(1/nYrs)
        # for making populations go extinct, we need to iterate across a function nYr times
        # each time we pass in the extant sites and get back the ones that survived to the next year
        #annual_survival <- function(occ, rate) occ[rbinom(n=length(occ),1,rate)==1]
        #occ <- annual_survival(occ, ann_persit_rate) # I can't get this to work inside sapply, but does work inside afor loop
        records <- list()
        for(i in 1:nYrs) {  #can start loop at zero
            for(j in 1:length(which.decline)){
                occ <- which(true_data[,which.decline[j]]==1) #site numbers where focal is present
                extinctions <- occ[rbinom(n=length(occ),1,ann_persit_rate[j])==0] # the site numbers at which extinctions will happen this year
                true_data[extinctions,which.decline[j]] <- 0 # set them to extinct
            }
            records[[i]] <- recording_cycle(pSVS=pSVS[1]+i*pSVS[2], true_data=true_data, max_vis=mv, VisRichSites=vrs, stochastic=stoch) #if starting loop at zero, use records[[i+1]]
        }
    }
    #records <- lapply(1:nYrs, FUN=recording_cycle, nVisits=nVisits, true_data=true_data) #simple version with constant number of visits
    
    records <- melt(records, id.vars=1:3) #simply appends the two list elements into a simple 
    names(records)[4] <- 'Year'
    return(records)
}


get_all_stats <- function(output, save_to_txt=T, sf=3){
    # 1/5/13: minor changes so that parameters are read from attr(output, ...) rather than passed
    
    if(dim(output)[3] == 1) {
        # 17/1/13: a quick hack to make the function work with only one replicate
        require(abind)
       
        # Added 9/5/3: doing this removes the attr
        simpars <- attr(output, "simpars")
        combos <- attr(output, "combos")
        id <- attr(output, 'id')
        time <- attr(output, "elapsed_time")
        
        output <- abind(output[,,1],output[,,1], along=3)

        # now add back the simpars
        simpars -> attr(output, "simpars")
        combos -> attr(output, "combos")
        id -> attr(output, 'id')
        time -> attr(output, "elapsed_time")

    }
    
    error_rates <- get_error_rates(output, save=F)
    summ_stats <- get_summary_stats(output, save=F, sf=sf)

    # how many replicates are valid - for Mixed model it's sometimes less than the specified
    p_val_cols <- grep('_p$',dimnames(output)[[1]]) #changed 16 Oct so only matches to end of string
    n_valid_reps <- apply(output[p_val_cols,,],1:2,function(a) sum(!is.na(a)))
    dimnames(n_valid_reps)[[1]] <- gsub(x=dimnames(n_valid_reps)[[1]], pattern='_p', repl='_nReps')
    
    x <- t(cbind(error_rates, summ_stats, t(n_valid_reps)))
    if(save_to_txt)	{
        id <- attr(output, "id")
        write.csv(x, file=paste('SimStats_',id,'_',dim(output)[3],'.csv', sep=''))
	} else {return(x)}
}

error_rate <- function(p, alpha=0.05) {p <- p[!is.na(p)]; length(p[p<=alpha])/length(p)}

test_power <- function(x, alpha=0.05, positive) {
    #modified from error_rate
    #now, the power is the proportion of replicates that are *both* less than alpha AND in the correct direction
    right_direction <- (x[1,] > 0) == positive #logical vector 
    p <- x[2,]
    if(any(p <= alpha & !right_direction, na.rm=T)) print('significant correlation opposite to true!')
    return(length(p[p<=alpha & right_direction])/length(p))
}


get_error_rates <- function(output, datecode=NULL, save_to_txt=T, sf=4, true_change=NULL){
	# first get the p-values from the methods that actually return p-values
	# 24 October: replaced dnorm with  one_to_two_tail(pnorm())): the former is conservative
	# 5-6 November: modified so that testing for power only looks for negative (or positive) associations (true_change)
	# 4/2/13 amended for cases where p includes NAs (typically MM under scenario C)
    # 1/5/13: minor changes so that parameters are read from attr(output, ...) rather than passed
    
    # is it power or validity? If not passed then read from output directly
    # fixed an error: 9 May so that power calculated in correct direction
    if(is.null(true_change)) true_change <- - attr(output, 'simpars')$decline

    #which columns contain the p-values
	p_val_cols <- grep('_p$',dimnames(output)[[1]]) #changed 16 Oct so only matches to end of string

	# if testing power, we need to know which columns match the p-value columns
	trend_names <- gsub('p$','trend',dimnames(output)[[1]][p_val_cols])
	trend_cols <- sapply(trend_names, grep, dimnames(output)[[1]])

	if(true_change==0) { 
		error_rates <- sapply(p_val_cols, function(i) apply(output[i,,], 1, error_rate)) #test validity
		plr = apply(one_to_two_tail(pnorm(output[dimnames(output)[[1]] == 'plr',,])),1,error_rate)
		Telfer = apply(one_to_two_tail(pnorm(output[dimnames(output)[[1]] == 'Telfer',,])),1,error_rate)
	} else { 
		error_rates <- sapply(1:length(p_val_cols), function(i) 
        apply(output[c(trend_cols[i], p_val_cols[i]),,], 2, test_power, positive=true_change>0)) #power (fixed error 21 Nov)
        plr = apply(one_to_two_tail(pnorm(output[dimnames(output)[[1]] == 'plr',,])),1,error_rate) #TO FIX - doesn't account for wrong direction
		Telfer = apply(one_to_two_tail(pnorm(output[dimnames(output)[[1]] == 'Telfer',,])),1,error_rate) #TO FIX - doesn't account for wrong direction
	}
	
	varnames <- dimnames(output)[[1]][p_val_cols]
	dimnames(error_rates)[[2]] <- substr(varnames, 1, regexpr('_',varnames) - 1)

	# now append the p-values from the standardised residuals
	error_rates <- cbind(prop.diff = NA, plr = plr, Telfer = Telfer,	error_rates)
    
    # for some combos (especially mixed model under scenario C) there are some reps which don't return enough data
    error_rates <- round(error_rates,sf)
    
    if (save_to_txt) write.csv(error_rates, file=paste('ErrorRates_',datecode,'_',dim(output)[3],'.csv', sep=''))
	return(error_rates)
	
	}


get_summary_stats <- function(output, datecode=NULL, sf=3, save_to_txt=T){
    # the various 'rows' of output contain a range of information types. Separate them out
    i <- grep('_trend', dimnames(output)[[1]]) # added 25/1/13
    i <- c(grep('Recs_', dimnames(output)[[1]]), i) # added 14/12/12, modified 22/1/13
    i <- c(grep('_t[0-9]', dimnames(output)[[1]]), i)
    i <- c(grep('pCombosUsed', dimnames(output)[[1]]), i) # if no mixed model this appends -1 (i.e. it gets ignored). Slight modification for multiple MMs (3/11/12)
    i <- c(grep('Fr_', dimnames(output)[[1]]), i) # if no Frescalo this appends -1 (i.e. it gets ignored). Added 17/1/13, modified 24/1

    # mean value for all stats (na.rm added 4/2/13)
    # 1/5/13: changed to median since Maes' method occasionally returns Inf (and is bounded at zero)
    x <- apply(output, 1:2, median, na.rm=T)[i,] 
    
    #strsplit(dimnames(x)[[2]], '_') # add this in future
    #x <- t(signif(x,sf)) # transpose and remove extraneous decimal places
    x <- t(round(x,sf)) # changed from above to avoid rounding the number of records to nearest 100!
    if (save_to_txt) write.csv(x, file=paste('SummaryStats_',datecode,'_',dim(output)[3],'.csv', sep=''))
    return(x)
}


incidental_records <- function(true_data, nYrs, nVisits, p_short, focal=F){
  # generates incidental records for just one species; either the focal or sppSVS
  # each year, generate a set of incidental records, all of the same species
  # first we need a vector of sites in which the target species is found
  # This doesn't work under test of power with nonfocal=T
  
  col <- 1+as.numeric(!focal) # column 1 for the focal, column 2 for sppSVS (the nonfocal)
  sites <- which(true_data[,col]==1) # set of sites where the target species if found
  
  if(length(p_short)==1){
    nInc <- p_short * nVisits # number of incidental records per year
    recs <- replicate(nYrs, sample(as.numeric(sites), nInc, repl=T)) #as.numeric is helpful to avoid confusion
    recs <- melt(recs)
    names(recs) <- c('Visit','Year','Site')
  } else {
    increment <- diff(unlist(p_short)) / nYrs # which extra ones should be included in each year
    nInc <- round((p_short$init + increment * 1:nYrs) * nVisits) # number of incidental records per year
    recs <- lapply(1:nYrs, function(i) sample(sites, nInc[i], repl=T))
    recs <- melt(recs)
    names(recs) <- c('Site','Year')
    recs$Visit <- melt(sapply(nInc, function(x) 1:x))$value
  }

  recs$Visit <- recs$Visit + nVisits # these are extra visits, over and above the normal
  recs$Species <- ifelse(focal, 'focal', 'sppSVS')
  return(recs[,c('Species','Visit','Site','Year')]) #structure matches main records, to which this appended
  }


iterate_all_scenarios <- function(nreps=1, nSites=1000, nSpecies=50, nYrs=10, pSVS=0.05, save_data=F, pFocal=list(Occ=0.5, DetP=0.5), vrs=F, mv=20, stoch=T,
                                  p_short=list(init=0,final=0.2), pDetMod=0.2, decline=0, id='', combos=F, Scenarios='BCDF', inclMM=F, Frescalo=F) {
    # this function is the wrapper for the standard procedure of generating data then analysing it and spitting the output
    # TO DO: allow save data to be an integer defining the number of reps to be printed. e.g. the first 100 out of 1000
    # As a workaround, I've added the 'inclMM' argument: if FALSE it bypasses the mixed model, which takes most time.
    # 8 April: I separated the data generation from data analsysis
    
    x <- system.time(
         output <- replicate(nreps, { # probeme when nreps>1. Hang on - is it when Frescalo=F?
            # first generate the records
            records <- generate_all_scenarios(nSites=nSites, nSpecies=nSpecies, nYrs=nYrs, pSVS=pSVS, save_data=save_data, 
                                              pFocal=pFocal, p_short=p_short, pDetMod=pDetMod, decline=decline, id=id, 
                                              combos=combos, Scenarios=Scenarios, mv=mv, vrs=vrs, stoch=stoch)

            # now analyse the data
            sapply(records, run_all_methods, inclMM=inclMM, Frescalo=Frescalo)
        })
    )
    attr(output,"simpars") <- list(nreps=nreps,nSpecies=nSpecies,nSites=nSites,nYrs=nYrs,pSVS=pSVS,pFocal=pFocal,p_short=p_short,
                                   pDetMod=pDetMod,decline=decline, visit_rich_sites=vrs, stochastic=stoch, max_vis=mv)
    attr(output,"combos") <- list(Scenarios=Scenarios, combos=combos, MM=inclMM, Frescalo=Frescalo)
    attr(output, "id") <- id
    attr(output, "elapsed_time") <- as.numeric(x[3])
    save(output, file=paste('SimOutput_',id,'_',dim(output)[3],'.rData', sep=''))
    return(output)
}

LenUniq <- function(x) length(unique(x)) #for calculating the list length from a list of species (allows duplicates, which are ignored)

median_occ <- function(true_data){
    # returns the column number of true_data corresponding to the species with the median distributions size (excluding the focal)
    # there's no which.median() function to this effectively creates that functionality
    x <- colSums(true_data)[-1]
    #x[order(x)] #lists species in increasing order
    name <- names(x[order(x)][ceiling(length(x)/2)])
    return(1 + as.numeric(gsub('spp','',name))) # add 1 because of the focal species
}


one_to_two_tail <- function(p) ifelse(p<0.5, p*2, (1-p)*2)


plot_error_rates <- function(output, error_rates, datecode){
	# the various 'rows' of output contain a range of information types. Separate them out
	p_val_cols <- grep('_p',dimnames(output)[[1]])
	recording_rate_cols <- grep('_t[1-9]',dimnames(output)[[1]])
	method_cols <- setdiff(1:dim(output)[1], c(p_val_cols, recording_rate_cols))
	
	# plot a histogram of the other parameters
	pdf(file=paste('Sim_T1Errors_',datecode,'_',dim(output)[3],'.pdf', sep=''), onefile=T, paper='a4', height=11.69, width=8.27)
	par(mfrow=c(dim(output)[2],length(method_cols)), mar=c(5,2,3,1), omi=c(0,0,0,0))

	for(j in 1:dim(output)[2]) # for each recording scenario
		for(i in 1:length(method_cols)) {# for each test statistic, excluding the p-values
			hist(output[method_cols[i],j,], xlab=dimnames(output)[[1]][method_cols[i]], col.lab='blue', main=NULL)
			if(i==1) {title(main=paste(dimnames(output)[[2]][j], sep=''), col.main='red')
			} else {title(sub=paste(quote(alpha),'=',error_rates[j,i], sep=''),col.sub='blue')}
			#} else {title(sub=paste(expression(alpha),'=',error_rates[j,i], sep=''),col.sub='blue')}
			#} else {title(sub=paste(substitute(alpha),'=',error_rates[j,i], sep=''),col.sub='blue')}
			#} else {title(sub=substitute(paste(alpha,'=',error_rates[j,i], sep='')),col.sub='blue')}
			#} else {title(sub=paste(quote(substitute(alpha)),'=',error_rates[j,i], sep=''),col.sub='blue')}
			#} else {title(sub=expression(paste(alpha,'=',eval(error_rates[j,i]))))}
			}	
	dev.off()	
	}


plot_error_rates_trace <- function(output, datecode=NULL, ignore.first=200){
	#plots how the error rates change with increasing replicate number
	#require(lattice)
	require(ggplot2)
    nreps <- dim(output)[3]
	if(nreps < ignore.first) sample_reps <- 2:nreps else sample_reps <- c(1:4*50, 250:nreps)
	er <- sapply(sample_reps, function(x) get_error_rates(output[,,1:x], save_to_txt=F, tr=attr(output, 'simpars')$decline), simplify=F)
	er <- melt(er)
	er$repnum <- sample_reps[er$L1]
	#pdf(file=paste('Sim_TraceErrors_',datecode,'_',dim(output)[3],'.pdf', sep=''), onefile=T, paper='a4r', width=11.69, height=8.27)
	#xyplot(data=er, value ~ repnum|Var1, gr=Var2, type='l', xlab='number of replicates', ylab='error rate', auto.key=T, scales=list(relation='free'), as.table=T)	
    names(er)[2] <- 'Method'
    qp <- qplot(data=er, x=repnum, y=value, col=Method, group=Method, geom='line') 
    qp + facet_wrap(~ Var1, nr=2) + xlab('number of replicates') + ylab('error rate')
    #dev.off() #pdf doesn't work!
	}


plot_summary_stats <- function(output, datecode){
	# the various 'rows' of output contain a range of information types. Separate them out
	# NB This will look congested if multiple mixed model thresholds were used
	t1_cols <- grep('_t1',dimnames(output)[[1]])
	t2_cols <- grep('_t2',dimnames(output)[[1]])
	#nc <- dim(output)[1]
	sc <- grep('MM_pCombosUsed', dimnames(output)[[1]]) # 'starting column' is given by this, the first summary statistic

	# plot a histogram of the other parameters
	pdf(file=paste('Sim_Summary_',datecode,'_',dim(output)[3],'.pdf', sep=''), onefile=T, paper='a4', height=11.69, width=8.27)
	par(mfrow=c(dim(output)[2],1+length(t1_cols)), mar=c(5,2,3,1), omi=c(0,0,0,0))

	for(j in 1:dim(output)[2]) # for each recording scenario
		for(i in c(sc,t1_cols)) {# for each summary stat
			hist(output[i,j,], xlab=dimnames(output)[[1]][i], col.lab='blue', main=NULL)
			if(i==sc) {title(main=paste(dimnames(output)[[2]][j], sep=''), col.main='red')}
			title(sub=paste('Mean=',round(mean(output[i,j,]),3), sep=''),col.sub='blue')
				}

	for(j in 1:dim(output)[2]) # for each recording scenario
		for(i in c(sc,t2_cols)) {# for each summary stat
			hist(output[i,j,], xlab=dimnames(output)[[1]][i], col.lab='blue', main=NULL)
			if(i==sc) {title(main=paste(dimnames(output)[[2]][j], sep=''), col.main='red')}
			title(sub=paste('Mean=',round(mean(output[i,j,]),3), sep=''),col.sub='blue')
				}
				
	dev.off()	
	}
	

plot_summary_stats3 <- function(output, datecode=NULL, maxN=7){
	# New version of above using boxplot over hist AND transposed
	# the various 'rows' of output contain a range of information types. Separate them out
	t1_cols <- grep('_t1',dimnames(output)[[1]])
	#t2_cols <- grep('_t2',dimnames(output)[[1]]) #should be just t1_cols + 1
	#nc <- dim(output)[1]
	#sc <- grep('MM_pCombosUsed', dimnames(output)[[1]]) # 'starting column' is given by this, the first summary statistic
	stat_names <- gsub('_t1', '', dimnames(output)[[1]][t1_cols])
	
	# plot a histogram of the other parameters
	pdf(file=paste('Sim_Summary_',datecode,'_',dim(output)[3],'.pdf', sep=''), onefile=T, paper='a4r', width=11.69, height=8.27)
	par(mfrow=c(length(t1_cols),7), mar=c(4,2,2,1), omi=c(0,0,0,0)) # 7 is the number of basic scenarios - any others get put on a second page
	
	for(i in 1:length(t1_cols)) # for each summary stat
		for(j in 1:maxN) {# for each recording scenario
			#temp <- output[t1_cols[c(i,i+1)],j,]
			temp <- output[c(t1_cols[i],t1_cols[i]+1),j,]
			dimnames(temp)[[1]] <- c('t1', 't2')
			boxplot(t(temp), xlab=stat_names[i], col.lab='blue', main=NULL)
			if(i==1) {title(main=paste(dimnames(output)[[2]][j], sep=''), col.main='red')}
			#title(sub=paste('Mean=',round(mean(output[i,j,]),3), sep=''),col.sub='blue')
			}

	dev.off()	
	}
	

recording_cycle_old <- function(nVisits, true_data, with_repl=T, max_vis=10) {
	# runs a 'recording cycle': nVisits are apportioned randomly among the sites (rows) in true_data: some sites get multiple visits
    # 3 December: modified to give 'fatter tails' in the distribution of recording effort among sites
    #   in UK recording data, about 1/4 of monads visited each year get further visits
    #   For the Dutch Dragonfly & Butterfly data, the figure is about 1/2
    #   The parameterization here is for UK.
    #   this makes the total number of visits and expectation, rather than fixed
    #  8/4/13: I removed the correlation between p(visited) and species richness
    # 24/4/13: New formulation 
    
    
    nS = nrow(true_data) # number of sites
  
    # this is the old version
    #sites_visited_old <- sample(1:nS, nVisits, repl=with_repl) #sites_visited is a vector of site identities
	# about 12% of sites receive multiple visits per year (for 250 visits to 1000 sites)

    # getting the number of visits per site from the binomial distribution: same distribution as above
    #vis_per_site <-rbinom(nS, max_vis, nVisits/(nS*max_vis)) # only 10% of visited sites get >1 visit
    #getting the number of sites from the betabinomial 
    #shape = 3000/nVisits # just about right to get 25% of visited sites having >1 visit
    #pr <- rbeta(nS,1,shape) # thanks to Steve Freeman for this insight
    #vis_per_site <- rbinom(nS, max_vis, pr*nVisits/nS) # increases to 25% the proportion visited twice or more
    #sites_visited_bb <- unlist(sapply(which(vis_per_site>0), function(i) rep(i, vis_per_site[i])))
    #problem - shape parameter means that nVisit is nonlinearly related to the acutal number of visits
    
    #14 December: revert to the old way, but make the probability of being selected related to richness 
    #rich <- rowSums(true_data)
    #temp <- rich - rpois(nS,2)
    #ss <- (temp/max(temp))^2 # makes negative values postive & increases the spread
    #sites_visited <- sample(1:nS, nVisits, repl=with_repl, prob=ss/sum(ss)) #sites_visited is a vector of site identities
    
    # 8 April: we no longer need 'memory' in the system
    sites_visited <- sample(1:nS, nVisits, repl=with_repl)
    
    #temp <- rich + rpois(nS,mean(rich)) - rpois(nS,mean(rich))
    #temp[temp<0] <- 0
    #sites_visited <- sample(1:nS, nVisits, repl=with_repl, prob=temp) #sites_visited is a vector of site identities
    
    # compare the differnet ways
    #stats <- list(sites_visited_old, sites_visited_bb, sites_visited)
    #sapply(stats, length) # number of visits
    #lapply(stats, function(x) table(table(x)))
    #sapply(stats, function(x) {x<- table(table(x)); sum(x)}) # number of sites
    #sapply(stats, function(x) {x<- table(table(x)); sum(x[-1])/sum(x)}) # prop getting multi visits
    #table(table(sites_visited))
    #the new way is better, but only 50% higher
	
    x <- visit_these_sites(sites_visited, true_data)
    records <- melt(x) #Var2 is the VisitNum
	names(records)[1:2] <- c('Species', 'Visit')
	records <- subset(records, subset=value==1)
	records$Site <- sites_visited[records$Visit]
	return(subset(records, select=c(Species, Visit, Site)))
}

recording_cycle <- function(pSVS=0.05, true_data, max_vis=10, VisRichSites=F, stochastic=F) {
    # runs a 'recording cycle': nVisits are apportioned randomly among the sites (rows) in true_data: some sites get multiple visits
    # 3 December: modified to give 'fatter tails' in the distribution of recording effort among sites
    #   in UK recording data, about 1/4 of monads visited each year get further visits
    #   For the Dutch Dragonfly & Butterfly data, the figure is about 1/2
    #   The parameterization here is for UK.
    #   this makes the total number of visits and expectation, rather than fixed
    #  8/4/13: I removed the correlation between p(visited) and species richness
    # 24/4/13: New formulation based on power law. NB this means number of visits varies from year to year
    #  2/5/13: I added back the correlatin between p(visited) and species richness, with an option to control
    #   I also added in the option for stochastic (the detault) vs determininistic numbers of visits    
    #  3/5/13: VisRichSites now additionally determines the number of visits received by each site (if selected)
    
    fx <- function(a=-3,b=-2, nSite=1000, maxV=10, stochastic=T){
        # This function defines the distribution of visits among sites in each year, based on the observed power law
        # the power law is truncated at maxV visits per year.
        
        # first define the log probabilities based on the model fit
        pr <- exp(a + b * log(1:maxV))
        
        # now append the probability of getting zero visits
        pr <- c(1-sum(pr), pr)
        
        #enS <- nSite * exp(y) ## the expected number of sites in each group
        # the number of visits at each site is defined by a multinomial distribution
        # this includes the number of sites with zero visits, so remove that from the calculation
        if(stochastic){
            n <- rmultinom(n=1, size=nSite, prob=pr)[-1]
        } else {
            n <- round(pr*nSite)
            # are there rounding errors?
            rem <- sum(n) - 1000 # corrected typo 18 june ('-' replaced '<')
            #if there are rounding errors, alter the number of sites in the zero category
            if(rem > 0) n[1] <- n[1] - rem
        }
        return(n)
    }
    
    fxp <- function(pSVS=0.05,prch=0.25, nSite=1000, maxV=10, stochastic=T){
        # a wrapper for fx in which the input parameters are on the measurement, rather than modelled scale:
        #   the proportion of sites receiving 1 visit & prch: the multiplier for doubling the number of visits 
        fx(a=log(pSVS), b=log2(prch), nSite=nSite, maxV=maxV, stochastic=stochastic)
    }
    
    # number of sites
    nS = nrow(true_data)
    
    # 3/5/13 VisRichSites gets used twice: if we only get one then use it twice
    if(length(VisRichSites)==1) {VisRichSites <- rep(VisRichSites, 2)
    } else {VisRichSites <- unlist(VisRichSites)}
    
    # get the number of sites receiving 1, 2, 3 visits
    #sv <- fxp(pSVS=pSVS, prch=0.25, nSite=nS, maxV=max_vis)
    sv <- fxp(pSVS=pSVS, prch=0.25, nSite=nS, maxV=max_vis, stochastic=stochastic) #corrected 18/6/13
    
    # set the probability that sites will be visited
    if(VisRichSites[1]) {prVis <- attr(true_data, "richness")
    } else {prVis <- NULL}

    # to figure out which sites get visited, sample sum(sv) sites without replacement
    # unlike sites_visited, each site appears only once in sites_to_visit 
    sites_to_visit <- sample(1:nS, sum(sv), repl=F, prob=prVis)
    
    # should we order the sites so that speciose sites ALWAYS get most visits?
    if(VisRichSites[2]) {
        # create a temporary dataframe of the selected sites and their richness
        nsp=attr(true_data, "richness")[sites_to_visit]
        sites_to_visit <- sites_to_visit[order(nsp)]
    }
    
    # unpack sv so that it has the same length as sites_to_visit
    nVis <- unlist(sapply(1:length(sv), function(i) rep.int(x=i, times=sv[i])))
    
    #visit_site_Ntimes(i=sites_to_visit[646], times=nVis[646], true_data)
    x <- mapply(FUN=visit_site_Ntimes, i=sites_to_visit, times=nVis, MoreArgs=list(true_data=true_data))
        
    records <- melt(x) #Var2 is the VisitNum - in the new version it goes from 1:max(nVis)
    names(records)[1:2] <- c('Species', 'RepVisit')
    records$Site <- sites_to_visit[records$L1]
    
    # renumber the visits so that they're unique within years, not within site:year combos
    uniq_vis <- dcast(records, Site + RepVisit ~ ., length)
    uniq_vis$id <- 1:nrow(uniq_vis)
    records <- merge(records, uniq_vis)
    records$Visit <- records$id
    
    # return the informative columns, without the occasions where a species was not observed 
    return(subset(records, subset=value==1, select=c(Species, Visit, Site)))
}

recording_visit <- function(spp_vector, p_obs=0.5, S=1){
	# simulate a visit to the site, OR a search for one species across multiple sites
	# takes a vector of species occurrences
	# the simplest way of doing this is to generate the observations for all sites first, even those without actual presences
    # in the default version, all species have equal probability
    # but species-specific detectabilities can be used also
    # detection probabilities can also be modified by the D parameter, which is the intensity of the search
    # for a complete search, S=1. (from Mark Hill)
    
    Hill_search <- function(pr, S){
        #takes a vector of species probabilities and the search intensity, S
        1 - exp(S * log(1-pr))
    }
    
    if (S==1) detected <- rbinom(n=length(spp_vector), size=1, prob=p_obs)
    else detected <- rbinom(n=length(spp_vector), size=1, prob=Hill_search(p_obs,S))
	# the realised set of observations is then the minimim value (i.e. a species was only recorded if both present and detected)
	apply(cbind(detected, spp_vector), 1, min)
}

resample <- function(x, ...) x[sample.int(length(x), ...)] # added 3/11 7 moved to separate function 5/2/13


run_all_methods <- function(records, min_sq=5, summarize=T, inclMM=2, Frescalo=TRUE){
	# 3 November: inclMM modified from a Boolean to the option of a number or vector of numbers, each of which specify the values of nsp for including
	# 6 November: added a 'stupid' model based on simply the number of visits & sites
	# 4 December: added 'pRecsBenchmark': the proportion of records made up by commonest 27%
    #takes a simulated dataset and runs each method
	
	######## two timeperiod models (0.01 seconds)
	splityr <- mean(range(records$Year))
	gridcell_counts <- Convert_records_to_2tp(records, splityr)	#get the number of sites in each time periods
	prop.diff <- with(subset(gridcell_counts, n1>=min_sq & n2 >0), (n2-n1)/n1)
	#plr <- rstandard(lm(log(n2) ~ log(n1), gridcell_counts, subset=gridcell_counts$n1 >= min_sq & gridcell_counts$n2 >0))
	plr <- rstandard(lm(log(n2) ~ log(n1), gridcell_counts, subset=n1 >= min_sq & n2 >0))
	Telfer <- fit_Telfer(gridcell_counts, min_sq)
	
    #we actually only want information from the focal species
	output <- cbind(prop.diff, plr, Telfer)[1,]
    
	# now fit the method of Maes et al 2012
    Maes <- fit_Maes(records, splityr, min_sq)[1,]
    output <- c(output, Maes_trend=Maes$trend, Maes_p=Maes$pval)
    
	attr(output, 'nSitesFocal') <- gridcell_counts[1,1:2]
	
    ######## DATA FRAME FOR visit-based analysis
	simdata <- cast_recs(records, resolution='visit') # 0.02 seconds
    
	######## WHICH SITES HAS THE FOCAL SPECIES EVER BEEN RECORDED ON? (added 1/5/13)
	focalsites <- unique(subset(records, Species=='focal')$Site)
        
	######## List Length (0.01 seconds)
	# 4/2/13: looking for an occasional bug
    x <- try(
        LL_model <- summary(glm(focal ~ Year + log2(L), binomial, data=simdata, subset = L>0))$coef
        , silent=T)
	if(class(x)=='try-error') save(records, file='LL1_try_error.rData')
	# end bug check
	output <- c(output, LL_trend=LL_model[2,1], LL_p=LL_model[2,4])
    
    ## a second version of the List Length in which only sites where the focal has been EVER recorded
	# added 1/5/13
    x <- try(
	    LL_model2 <- summary(glm(focal ~ Year + log2(L), binomial, data=simdata, 
                                subset = L>0 & Site %in% focalsites))$coef
	    , silent=T)
	if(class(x)=='try-error') save(records, file='LLfs_try_error.rData')
	# end bug check
	output <- c(output, LLfs_trend=LL_model2[2,1], LLfs_p=LL_model2[2,4])
    
	######## Ball method 
	# modified 1/5/13 to include multiple versions of this method
    
    # first, the simple version (as used previously)
	num_visits_recorded <- with(simdata, tapply(focal, Year, sum))
	num_visits_that_year <- with(simdata, tapply(Site, Year, length))
	failures <- num_visits_that_year - num_visits_recorded#[,2] #col 2 is focal species (col 1 is year)
	Year <-unique(simdata$Year)
		
	x <- try(# 4/2/13: looking for an occasional bug
	    Ball_model <- summary(glm(cbind(num_visits_recorded,failures) ~ Year, family=binomial))$coef
	    , silent=T)
	if(class(x)=='try-error') save(records, file='try_error.rData')
	# end bug check
	output <- c(output, VRsimple_trend=Ball_model[2,1], VRsimple_p=Ball_model[2,4])    
    
    # Note that in the hoverfly atlas, Ball fits the logit(proportion), with number of hectads & visits as covariates
    # the binomial version of this model should render these covariates unnecessary
    # he also includes number of unique recorders
    
    # Next, we use a version in which only sites where the focal was ever recorded
    # this is a filter discussed by Stuart Ball but not described in the hoverfly atlas
    # Ball calls these visits 'suitable', although they would be defined at a coarser resolution than used here.
    #num_visits_recorded <- aggregate(simdata$focal, by=list(simdata$Year), FUN=sum)
	num_visits_recorded <- with(subset(simdata, Site %in% focalsites), tapply(focal, Year, sum))
	num_visits_that_year <- with(subset(simdata, Site %in% focalsites), tapply(Site, Year, length))
	failures <- num_visits_that_year - num_visits_recorded#[,2] #col 2 is focal species (col 1 is year)
	YearFS <-unique(subset(simdata, Site %in% focalsites)$Year) # modified 9/5 to catch years when no focal sites are visited
    
	x <- try(# 4/2/13: looking for an occasional bug
	    Ball_model <- summary(glm(cbind(num_visits_recorded,failures) ~ YearFS, family=binomial))$coef
	    , silent=T)
	if(class(x)=='try-error') save(records, file='VR_try_error.rData')
	# end bug check
    output <- c(output, VRfs_trend=Ball_model[2,1], VRfs_p=Ball_model[2,4])

    # Finally, we'll use the version in the hoverfly atlas
    # this is a logistic regression with covariates. 
    # Stuart Ball told me (8/5/13) that he did this just because it was easier to explain
    # but it's not clear to me whether we should restrict the data to Years when focal sites were visited
    nRecsFocal <- with(subset(records, Species=='focal'), tapply(Species,Year,length))
    recs_temp <- subset(records, Year %in% YearFS)
    nRecsTotal <- with(recs_temp, table(Year))
	n_Sites <- with(recs_temp, tapply(Site, Year, LenUniq))
	n_Vis <- with(recs_temp, tapply(Visit, Year, LenUniq))

    # 14/5/13: there's an occasional problem here
    # in some datasets, there are zero focal records for a particular year
    # this means that the lengths of the vectors are different
	if(length(nRecsFocal) < length(nRecsTotal)){
	    nRecsFocal <- nRecsFocal[match(dimnames(nRecsTotal)[[1]], dimnames(nRecsFocal)[[1]])]
	    nRecsFocal[is.na(nRecsFocal)] <- 0.01 # a fudge to prevent NA/-Inf in the logit
	    }
    
    x <- try({
        pF <- nRecsFocal/nRecsTotal
	    logit_pF <- log(pF/(1-pF))
        Ball_model <- summary(lm(logit_pF ~ YearFS + n_Vis + n_Sites))$coef
	    }, silent=T)
    if(class(x)=='try-error') save(records, file='VRhov_try_error.rData')

    output <- c(output, VRhov_trend=Ball_model[2,1], VRhov_p=Ball_model[2,4])    
	
    
    ######## Mixed model: 0.16 seconds to cast the data, 0.3 to fit the model
	
    # 21 June: add the MM without any kind of threshold
    # this makes it the same as Visit rate but with a site effect
    # we're actually including it in order to verify that Occupancy will be ok under B3f
    MM <- fit_ladybirdMM_bin(simdata, nsp=1, nyr=1)[c(1,4,5)]
	names(MM) <- c('trend', 'p', 'pCombosUsed')
	names(MM) <- paste('MMbin0_', names(MM), sep='')
	output <- c(output, MM) #     
    
    #MMdata <- cast_recs(records, resolution='kmyr')
	#i <- MMdata$L >= 2
	#i[i==T] <- is.gridcell.wellsampled(MMdata$Site[i], n=3)
	#MM <- glmer(focal ~ Year + (1|Site), MMdata, subset=i, family=binomial)
	#MM <- as.numeric(summary(MM)@coefs[2,]) # keep this a two step process in case we later decide to extract other info
	#output <- c(output, MM_trend=MM[1], MM_p=MM[4], MM_pCombosUsed=sum(i)/length(i)) # M_pVisitsUsed added: 16 October
	if(sum(inclMM) > 0){ # option to bypass fitting the mixed model, e.g. if set to F
	    MMdata <- cast_recs(records, resolution='kmyr')
	    for (i in inclMM) {
		    # the standard MM, as year-monad (ym) resolution
            MM <- fit_LadybirdMM(MMdata, nsp=i, nyr=3)[c(1,4,5)]
			names(MM) <- c('trend', 'p', 'pCombosUsed')
			names(MM) <- paste('MMber', i,'sp_', names(MM), sep='')
			output <- c(output, MM) #
            
            # repeat using disaggregated (visit) data (Added 24/1/13)
            # NB Currently the nyr parameter actually filters on the number of visits, not unique years
            #MM <- fit_LadybirdMM(simdata, nsp=i, nyr=3)[c(1,4,5)]
            #names(MM) <- c('trend', 'p', 'pCombosUsed')
            #names(MM) <- paste('MMv', i,'sp_', names(MM), sep='')
            #output <- c(output, MM) #    
            
            # instead of disaggregated visit data, use the Binomial model
            # it should be identical to the MMv model, but with greater numerical stability & quicker
            # it's more robust than the bernoulli version above (previously MMym)
            MM <- fit_ladybirdMM_bin(simdata, nsp=i, nyr=3)[c(1,4,5)]
            names(MM) <- c('trend', 'p', 'pCombosUsed')
            names(MM) <- paste('MMbin', i,'sp_', names(MM), sep='')
            output <- c(output, MM) #    
	}}

    ######## Frescalo (by Tom August & Colin Harrower)
	if (sum(Frescalo)>0){
	 # Set directory where Frescalo is located
	 	  # Set path to Frescalo and output folder, which is platform dependent 
	  if (grepl("linux", R.version$platform)){
	    #frescalo_path <- '/users/hails/tomaug/NEC04273/BSBI Redlisting/Master Code/Frescalo/Frescalo_V3/Frescalo_3a_linux.exe'
	    #output_dir <- '/users/hails/tomaug/NEC04273/BSBI Redlisting/Master Code/Frescalo/Output'	
		frescalo_path <- '/prj/NEC04273/BSBI Redlisting/Master Code/Frescalo/Frescalo_V3/Frescalo_3a_linux.exe'
	    output_dir <- '/prj/NEC04273/BSBI Redlisting/Master Code/Frescalo/Output'
	  } else { 
	    frescalo_path <- "P:/NEC04273_SpeciesDistribution/Workfiles/Range change sims/Frescalo/Frescalo files/Frescalo_3.exe"
	    output_dir <- "P:/NEC04273_SpeciesDistribution/Workfiles/Range change sims/Frescalo/Output"
	  }
	 
	 for (i in Frescalo) {
    #i=1: 1 year per time period
    #i=2: 2 time periods (splitting years in two)
	  Tfac_Stdev <- sims_to_frescalo(records, output_dir=output_dir, frescalo_path=frescalo_path, tp2=as.logical(i-1))
    #get p-value & trend from TimeFactors
	  #frescalo_summ <- frescalo_bootstrap(Tfac_Stdev) #bootstrapping nearly ALWAYS returns a significant p-value
	  frescalo_summ <- frescalo_trend(Tfac_Stdev)
	  names(frescalo_summ) <- c('trend','p')
	  names(frescalo_summ) <- paste('Frescalo', i, 'tp_', names(frescalo_summ), sep='')
      output <- c(output, frescalo_summ) 
	}}
  
	######## Stupid method - a poisson regression of number of visits & sites
	NRecMod <- summary(glm(num_visits_recorded ~ as.numeric(dimnames(num_visits_recorded)[[1]]), 'poisson'))$coef
	output <- c(output, nRecords_trend=NRecMod[2,1], nRecords_p=NRecMod[2,4])
	
	num_sites_recorded <- with(subset(records, Species=='focal'), tapply(Site, Year, LenUniq))
	NSmod <- summary(glm(num_sites_recorded ~ as.numeric(dimnames(num_sites_recorded)[[1]]), 'poisson'))$coef
	output <- c(output, nSites_trend=NSmod[2,1], nSites_p=NSmod[2,4])
	
	#append summary info this to the output (doesn't work as attr when looping over multiple datasets
	if (summarize) {
        recording_summary <- summarize_recording(simdata, length(levels(records$Species))) # 0.03 seconds
       
        if (sum(Frescalo)>0){ # ~NB if(length(Frescalo) >1) then the following lines will only get info from the last run 
            if (grepl("linux", R.version$platform)){
                #frst <- read.csv("/users/hails/tomaug/NEC04273/BSBI Redlisting/Master Code/Frescalo/Frescalo_V3/Stats.csv")
			    frst <- read.csv("/prj/NEC04273/BSBI Redlisting/Master Code/Frescalo/Frescalo_V3/Stats.csv")			
			}else{
			    frst <- read.csv("P:/NEC04273_SpeciesDistribution/Workfiles/Range change sims/Frescalo/Frescalo files/Stats.csv")    
			}
        output <- c(output, Fr_Phi= attr(Tfac_Stdev, 'Phi'), Fr_MedianAlpha = median(frst$Alpha))
		}            
        
        #get information on the prevelance of benchmark species
        # NI 24/1/12
        # The stats below are based on numbers of records
        # In theory, it should be possible to extract the equivalent info from the Frescalo trend file
        # I'll wait on this until Tom August has fixed the column names of the different versions
        rps <- with(records, table(Species))
        bench_thresh <- quantile(rps, 0.73)
        pRecsBench <- sum(rps[rps >= bench_thresh])/sum(rps)
                
        output <- c(output, Recs_pBnchmk=pRecsBench, Recs_qFocal=mean(rps>rps[1]), recording_summary)
	}
    return(output)
}


sims_to_frescalo <- function (records,output_dir,frescalo_path,tp2=FALSE){
  # Tom August - 06/11/12
  # original is in 'Frescalo/Process files/sims_to_frescalo.r'
  # take columns Frescalo needs
  cols<-c("Site","Species","Year")
  fres_in<-records[cols]
  
  # bin years into two time periods
  if(tp2){
    split=mean(unique(fres_in$Year))
    lower<-fres_in$Year<=split # NI added '=' sign to ensure no data is lost
    upper<-fres_in$Year>split
    fres_in[lower,]$Year<-1
    fres_in[upper,]$Year<-2
  }
  fres_in<-unique(fres_in) # added 17/1/13 to remove duplicate entries
  # call the frecalo function
  Tfac_StDev<-run_fresc_sims(fres_in,output_dir,frescalo_path)
    
  # return the time factor and standard deviation for the species 'focal'
  return(Tfac_StDev)
}


save_records <- function(records, id='', as.R.obj=T) {
	#saves the output in one of two formats (csv works but has not been fully tested: two empty files in place of true_data)
	#ideally replace timestamp with a replicate number
    #25/4/13: Changed to include 'id', which by default also includes Datestamp
	#names(x)[[1]] <- 'TrueData'
	#names(x) <- gsub('_', '', names(x))
	path=paste('SimRecords', sep='')
	#datecode <- format(Sys.Date(),'%y%m%d')
	timestamp <- format(Sys.time(),'%H%M%S')
	if(as.R.obj) {save(records, file=paste(path,'/SimData_',id,'_',timestamp,'.rData', sep=''))
	} else {for (i in 1:length(x)) write.csv(x[[i]], file=paste(path,'/SimData_',timestamp,'_',names(x)[i],'.csv', sep=''))}
}


shorten_lists <- function(records, p_short=0.25){ #defaults to 25% of lists being short
	# a clever way to simulate short lists is to delete records from 'well-sampled' set by stratified random sampling
	# p_short defines the proportion of lists that should be short
	# If a single number then it's constant across years.
	# if numeric, this remains constant over time. If two numbers, the second number identifies the FINAL value
	increment <- ifelse(length(p_short)==1, 0, diff(unlist(p_short)) / max(records$Year)) # which extra ones should be included in each year
	
	#df of all the visits (Site is redundant but worth retaining for clarity)
	uniq_vis  <- dcast(records, Year + Site + Visit ~.,length, value.var='Species') #0.06 seconds
	uniq_vis$visit_id <- 1:nrow(uniq_vis) # this is necessary to make x easy to understand
	
	prop_short <- p_short[[1]] + increment * uniq_vis$Year
	uniq_vis$target_LL <- targetLL(prop_short)
	uniq_vis$recs_to_remove <- uniq_vis$'NA' - uniq_vis$target_LL # 'NA' is the produced by dcast(): it's the list length for each visit
	uniq_vis <- subset(uniq_vis, subset=recs_to_remove>0) #restrict uniq_vis to those where pruning is necessary

	records$record_id <- 1:nrow(records)
	records_tmp <- merge(records, uniq_vis, all=F) # records in visits to be 'pruned'

	#~0.1 secs to here: the last bit is quite slow! (0.1-1.4 seconds, depending on nVisits & p_short)
	#recs_to_omit <- unlist(apply(uniq_vis[,c('visit_id','recs_to_remove')], 1, function(x) { #for each visit that has been selected
	#	recs <- subset(records_tmp, visit_id==x[1])$record_id #the identities of all the records in selected visits
	#	sample(recs, size=x[2])
	#	})) #names(recs_to_omit) contains the visit number plus the number within the sample, ie. '23' is the 3rd visit to be sampled from the visit 2 
	
	# this version is much quicker, and *much* less dependent on p_short
	#recs_to_omit <- sapply(1:3, function(i) unlist(with(subset(records_tmp, recs_to_remove==i), tapply(record_id, visit_id, sample, i))))
	#26 October: corrected error in version above!
	
	recs_to_omit <- sapply(1:max(uniq_vis$recs_to_remove), function(i) with(subset(records_tmp, recs_to_remove==i), tapply(record_id, visit_id, resample, i)))
	
	# so the pruned data look like this:
	records <- subset(records, !record_id %in% unlist(recs_to_omit), select=-ncol(records)) # added bug fix 'unlist' 25/10
	#return(records)
	}


subset_recs_by_X <- function(records, X, init=0.5){
	# a clever way to simulate increasing numbers of visits OR Sites over time is to take the well-recorded set and strip out the data based on a set of rules
	# X should be either the vector of sites or visit identities
	# init defines the number or proportion to be surveyed in year zero
    # The subset is made easy because the visits and sites are numbered sequentially
    # sort added 25/1/13 for safety!
    # modified 25/4/13 to account for the fact that nVisits now varies among years

    # first calculate the number of visits in the each year
    nV <- with(records, tapply(Visit, Year, LenUniq))
    #nV <- with(records, tapply(Visit, Year, max)) # should be same as above, but larger by 1 in some years!
    
    #if(init < 1) init <- init * max_vis #max(records$Visit)
    increment <- (1 - init) / max(records$Year) # additional proportion of visits to be included in each year
    pV_include <- init + sort(unique(records$Year)) * increment
    nV_include <- round(nV * pV_include) # vector of length nYear defining the number of records to include
      
    i <- X <= nV_include[records$Year]
	#table(records$Year, i) # looks sensible
	return(records[i,])
}


subset_recs_Biased <- function(records, init=0.5, bias=2, target=NULL){
    # 25 January 2013
    # Modified version of subset_recs_by_X() in which sampling is biased in favour of sites where a target species is present
    # The parameter target is a vector of species presence-absence (either the focal or non-focal)
    # The default setting is for bias=2, indicating that sites with the target species are 2x likely to be visited
    # This parameter is the odds ratio for selecting sites where Target is present compared with where it's absent
    # this is actually over-engineered since the number of visits in the original (unsubsetted) records is fixed across years
    # 8/4/13: added a 'select' argument to the subset statement so that records have same strucutre as other Scenarios
    # modified 25/4/13 to account for the fact that nVisits now varies among years
    
    # first calculate the number of visits in the each year
    nV <- with(records, tapply(Visit, Year, LenUniq))
    #nV <- with(records, tapply(Visit, Year, max)) # should be same as above, but larger by 1 in some years!
    
    #if(init < 1) init <- init * max_vis #max(records$Visit)
    increment <- (1 - init) / max(records$Year) # additional proportion of visits to be included in each year
    pV_include <- init + sort(unique(records$Year)) * increment
    nV_include <- round(nV * pV_include) # vector of length nYear defining the number of records to include

    #df of all the visits
    uniq_vis  <- dcast(records, Year + Site + Visit ~.,length, value.var='Species') #0.06 seconds
    # the NA column is the List Length
    uniq_vis$visit_id <- 1:nrow(uniq_vis) # this is necessary to make x easy to understand
    
    # add the unique visit identity to records
    records <- merge(records, uniq_vis, by=c('Site','Year','Visit'))    
    
    # now create a vector of length = nrow(uniq_vis) indicating whether the focal species is present
    uniq_vis$target <- target[match(uniq_vis$Site, as.numeric(gsub('site','',names(target))))]

    # convert these into probabilities of being sampled, using bias as the odds ratio
    uniq_vis$p_sel <- ifelse(uniq_vis$target, bias, 1)
      
    vis_to_include <- unlist(sapply(unique(records$Year), 
        function(i) { 
            vis_tmp <- subset(uniq_vis,Year==i) # the subset of visits for that year
            # this sometimes fails.
            # this can happen when the number of actual visits is less than expected
            # in year 10 we're sampling the expected number of visits from a vector that should be the same length (e.g. 250)
            # occasionally it's slightly shorter. Possible reasons: one visit had zero species.
            # solution: if statement to bypass
            if(nrow(vis_tmp) <= nV_include[i]) x <- vis_tmp$visit_id
            else {  
                x <- try(sample(vis_tmp$visit_id, size=nV_include[i], replace=F, prob=vis_tmp$p_sel), silent=T)
                if(class(x)=='try-error') {
                    error_report <- list(records=records, uniq_vis=uniq_vis, include=nV_include, i=i)
                    save(error_report, file='try_error.rData')
                }}
            return(x)
        }))
        
    recs_to_return <- subset(records, visit_id %in% vis_to_include, select=c('Species','Visit','Site','Year'))
    # check the output
    #with(uniq_vis, plot(tapply(target, Year, sum), ylim=c(0,max(Visit)), type='l', col='red')) # even across years
    #with(subset(uniq_vis, visit_id %in% vis_to_include), lines(tapply(target, Year, sum), ylim=c(0,max(Visit)))) # increasing over time

    return(recs_to_return)
}


summarize_recording <- function(simdata, nsp){  
  ## Calculate the basic stats, each of the two time periods (split the years evenly)
  # get the number of nonfocal secies overall and in each visit
  
    simdata$nonfocal <- (simdata$L - simdata$focal)/(nsp-1) # proportion of species per visit
    splityr <- mean(range(simdata$Year))
  
  # stats on the visit rate (number of sites visited per year, and multiple times per year)
  Site_visit_stats <- acast(simdata, Site ~ Year, fun=LenUniq, value.var='Visit')	
  summary_stats <- c(
    Recs_tot_visits = nrow(simdata),
    Recs_Sp_Yr = sum(simdata$L)/(nsp * max(simdata$Year)), #added 9/5/13
    Recs_total = sum(simdata$L),
    Recs_focal = sum(simdata$focal),
    VisPerYr_t1 = mean(colSums(Site_visit_stats)[1:splityr]),
        #mean(with(simdata, tapply(Visit, Year, LenUniq))[1:splityr]), # changed now that 'Visit' is not unique
    VisPerYr_t2 = mean(colSums(Site_visit_stats)[-c(1:splityr)]),
        #mean(with(simdata, tapply(Visit, Year, LenUniq))[-c(1:splityr)]),
    focal_t1 = mean(subset(simdata, Year<=splityr)$focal),
    focal_t2 = mean(subset(simdata, Year>splityr)$focal),
    nonfocal_t1 = mean(subset(simdata, Year<=splityr)$nonfocal),
    nonfocal_t2 = mean(subset(simdata, Year>splityr)$nonfocal),
    meanL_t1 = mean(subset(simdata, Year<=splityr)$L),
    meanL_t2 = mean(subset(simdata, Year>splityr)$L),
    p_shortLists_t1 = mean(subset(simdata, Year<=splityr)$L <= 3), #proportion of lists that are 'short' (i.e. less than 4)
    p_shortLists_t2 = mean(subset(simdata, Year>splityr)$L <= 3),
    MeanSitesVisited_t1 = mean(colSums(Site_visit_stats>0)[1:splityr]),
    MeanSitesVisited_t2 = mean(colSums(Site_visit_stats>0)[-c(1:splityr)]),
    #MeanMultiVisitSites_t1 = mean(colSums(Site_visit_stats>1)[1:splityr]),
    #MeanMultiVisitSites_t2 = mean(colSums(Site_visit_stats>1)[-c(1:splityr)])
    PrSitesMultiVisit_t1 = mean((colSums(Site_visit_stats>1)/colSums(Site_visit_stats>0))[1:splityr]),
    PrSitesMultiVisit_t2 = mean((colSums(Site_visit_stats>1)/colSums(Site_visit_stats>0))[-c(1:splityr)])
    )
  return(summary_stats)	
}


targetLL <- function(p_short, ratio=c(4,2,1.5)) {
	# p_short defines the proportion of short lists 
	# balance  defines how that proportion is split between lists of 1, 2 and 3 species long
	# so we draw a number for the number of target species, with NA representing all lists longer than 3
	# this version will return a vector
    # 3 December: I changed the 'balance' argument to reflect better the balance of real recording scheme data
    balance <- ratio/sum(ratio)
	sapply(p_short, function(x) sample(c(1:3,NA), size=1, prob=c(balance*x, 1-x)))
	}
	

visit_site_Ntimes <- function(i, times, true_data){
    # 24 April: new wrapper for recording_visit
    replicate(times, recording_visit(true_data[i,], p_obs=attr(true_data, "p_detect")))
}

visit_these_sites <- function(sites_visited, true_data){
	# takes a vector of sites to be visited and the set of actual observations
	# 3 December: I added the 'pr' argument to allow different detection probabilities among species 
	sapply(sites_visited, function(i) recording_visit(true_data[i,], p_obs=attr(true_data, "p_detect")))
}

#end