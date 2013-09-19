# Nick Isaac 7 May 2013
# Extract and explore the simulation results
# this time we'll collate the output

# Updated 29 May to include Power & Occupancy results
# Updated 22 july to add the mixed model with no threshold for data
# updated 13 Sept with the corrected version of Maes' model

library(ggplot2)
library(reshape2)

homedir <- getwd()

source('extract_collate_stats.R')
source('Sim_functions.r')

datecode <- format(Sys.Date(), '%y%m%d')

###################### Load the results and point to an output folder

sinkdir <- 'results/results 130912'

# first load the production set of results
nrow(output0 <- extract_output('results/results 130509')) #1155000

# now load the new results: these should be on identical data but only using the mixed model with 0 threshold
nrow(output1 <- extract_output('results/results 130621 +MM0')) #882000

# the 130912 replaces 130621, with a corrected version of the Maes results
nrow(output2 <- extract_output('results/results 130912')) #882000

###################### 
# test the various sets of output
temp <- output1$value==output2$value
table(temp) # only 50% match
# for(i in 1:ncol(output1)) print(table(output1[,i]==output2[,i])) # all others TRUE except id, which is good
table(temp[-grep('Maes',output1$Statistic)]) # 370644 mismatches!
# which ones are they?
mismatches <- setdiff(which(temp==F), grep('Maes',output1$Statistic))

plot(output1$value[mismatches], output2$value[mismatches]) # very close to 1:1

hist(output1$value[mismatches] - output2$value[mismatches]) # ahh 
summary(output1$value[mismatches] - output2$value[mismatches]) # a couple of extremes, but otherwise rounding errors 

extremes <- abs(output1$value - output2$value) > 0.1
temp <- cbind(subset(output1[extremes,], !grepl('Maes',Statistic), select=1:4),
  subset(output2[extremes,], !grepl('Maes',Statistic), select=4))
dim(temp) # 25 cases
temp # all where the mixed model with no thresholds gives different p-values (out of 3500).
# so we can be confident that the rest of the results are the same (due to setting the randseed constant)

######################Got to here

# stick them together
output <- droplevels(rbind(
            subset(output0, grepl('_p$',Statistic) & !grepl('Maes', Statistic)), # exclude Maes
            subset(output2, grepl('_p$',Statistic) & grepl('MMbin0', Statistic)), # just the MM0 results
            subset(output2, grepl('_p$',Statistic) & grepl('Maes', Statistic)) # just the MM0 results
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
# the new occupancy results are in a set of six summary files (plus the full output in many files)

occfiles <- list.files('results/Occupancy 130730', pattern='NumberOfType')
lapply(occfiles, )


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


# for Occupancy, we want to know what happens just under B3
png(paste0(sinkdir,'/validityB3 ',datecode,'.png'), wi=500, hei=500)
qp <- qplot(data=subset(ER, decline==0 & grepl('B3f', Scenario) & Method %in% c('VRsimple', 'Occupancy', 'MMbin0', 'LL')),
      x=pSVS, y=ER, geom='line', col=Method, size=I(1))         
qp + geom_hline(yintercept =0.05)
dev.off()

######################################################################## NOW REPEAT FOR POWER

# calculate the Type 2 error rate
ER$T2 <- ifelse(ER$decline != 0, 1 - ER$ER, NA)

# build the ggplot object
ggP <- ggplot(data=subset(ER, !is.na(T2)), aes(x=pSVS, y=T2, gr=Method))
ggP <- ggP + geom_line(aes(col=Method), size=1)
#ggP <- ggP + scale_y_log10()
ggP <- ggP + ylab('Type II error rate') + xlab('recording intensity')
ggP <- ggP + facet_wrap(~Scenario, nrow=1)
ggP


# MAKE SOME PLOTS OF VALIDITY
png(paste0(sinkdir,'/PowerA ',datecode,'.png'), wi=500, hei=500)
ggP %+% subset(ER, grepl('A', Scenario)) # All scenario A
dev.off()

png(paste0(sinkdir,'/Power1 ',datecode,'.png'), wi=1024, hei=500)
ggP # All  
dev.off()

png(paste0(sinkdir,'/Power2 ',datecode,'.png'), wi=1024, hei=500)
ggP %+% subset(ER, Scenario != 'C1_pShortLEven'
               & !Method %in% c('nRecords', 'Maes', 'VRhov')) # all except the really bad methods
dev.off()

png(paste0(sinkdir,'/PowerMM ',datecode,'.png'), wi=1024, hei=500)
ggP %+% subset(ER, grepl('MM', Method))# just MM
dev.off()

png(paste0(sinkdir,'/Power3 ',datecode,'.png'), wi=1024, hei=500)
ggP %+% subset(ER, Scenario != 'C1_pShortLEven'
               & !Method %in% c('nRecords', 'Maes')
               & !grepl('VR', Method)
               & !grepl('LL', Method) 
               & !grepl('A', Scenario)
               & !grepl('MMber', Method)
)# the good methods
dev.off()




