# Nick Isaac 7 May 2013
# Extract and explore the simulation results
# this time we'll collate the output

# Updated 29 May to include Power & Occupancy results
# Updated 22 july to add the mixed model with no threshold for data

library(ggplot2)
library(reshape2)

homedir <- getwd()

source('extract_collate_stats.R')
source('Sim_functions.r')

datecode <- format(Sys.Date(), '%y%m%d')

######################

sinkdir <- 'results/results 130621 +MM0'

# first load the production set of results
nrow(output0 <- extract_output('results/results 130509')) #1155000

# now load the new results: these should be on identical data but only using the mixed model with 0 threshold
nrow(output1 <- extract_output('results/results 130621 +MM0')) #882000

output <- droplevels(rbind(
            subset(output0, grepl('_p$',Statistic)),
            subset(output1, grepl('_p$',Statistic) & grepl('MMbin0', Statistic))
            ))
output$Method <- gsub(x=output$Statistic, pattern='_p', repl='')

###################################################################################### ERROR RATES
# calculate the error rates
ER <- dcast(output, Method + Scenario + pSVS + decline ~ ., fun=error_rate, value.var='value')
names(ER)[ncol(ER)] <- 'ER' # default is to call it 'NA'

# get the sample sizes, just for a sense-check
ER$nreps <- as.numeric(acast(output, Method + Scenario + pSVS + decline ~ ., fun=length, value.var='value'))

###################################################################################### OCCUPANCY
# Now import the Occupancy results and match them up.
# the results file is modified into a readable format from the spreadsheet sent to me by Marnix
occstats <- read.csv('results/results 130509/Occupancy results 130527.csv')
occstats$Method <- 'Occupancy'
occstats$ER <- with(occstats, X..sig.trends/nreps)

# merge
#ER <- rbind(ER, subset(occstats, select=intersect(names(occstats), names(ER))))
#ER <- rbind(ER, subset(occstats, select=names(occstats) %in% names(ER)))
ER <- rbind(ER, occstats[,names(occstats) %in% names(ER)])

ER$intensity <- factor(ER$pSVS, levels=c(0.05,0.07,0.1), labels=c('low','med','high'), ordered=T)

dim(ER) # 630 7   

######################################################################## PLOTTING VALIDITY


# build the ggplot object
ggV <- ggplot(data=subset(ER, decline==0), aes(x=pSVS, y=ER, gr=Method))
#ggV <- ggplot(data=subset(ER, decline==0), aes(x=intensity, y=ER))#, gr=Method))
ggV <- ggV + geom_line(aes(col=Method), size=1)
ggV <- ggV + scale_y_log10() + geom_hline(yintercept =0.05)
ggV <- ggV + ylab('Type I error rate') + xlab('recording intensity')
#ggV <- ggV + geom_point(aes(col=Method, size=3))
ggV <- ggV + facet_wrap(~Scenario, nrow=1)
ggV


# MAKE SOME PLOTS OF VALIDITY
png(paste0(sinkdir,'/validityA ',datecode,'.png'), wi=500, hei=500)
ggV %+% subset(ER, decline==0 & grepl('A', Scenario)) # All scenario A
dev.off()

png(paste0(sinkdir,'/validity1 ',datecode,'.png'), wi=1024, hei=500)
ggV # All validity 
dev.off()

png(paste0(sinkdir,'/validity2 ',datecode,'.png'), wi=1024, hei=500)
ggV %+% subset(ER, decline==0 &  Scenario != 'C1_pShortLEven'
               & !Method %in% c('nRecords', 'Maes', 'VRhov')) # all except the really bad methods
dev.off()

png(paste0(sinkdir,'/validityMM ',datecode,'.png'), wi=1024, hei=500)
ggV %+% subset(ER, decline==0 & grepl('MM', Method))# just MM
dev.off()

png(paste0(sinkdir,'/validity3 ',datecode,'.png'), wi=1024, hei=500)
ggV %+% subset(ER, decline==0 & Scenario != 'C1_pShortLEven'
               & !Method %in% c('nRecords', 'Maes')
               & !grepl('VR', Method)
               & !grepl('LL', Method) 
               & !grepl('A', Scenario)
               & !grepl('MMber', Method)
               )# the good methods
dev.off()

############ NOW REPEAT FOR POWER
