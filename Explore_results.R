# Nick Isaac 7 May 2013
# Extract and explore the simulation results

library(ggplot2)
library(reshape2)

homedir <- getwd()

source('extract_collate_stats.R')

datecode <- format(Sys.Date(), '%y%m%d')

###################### Load the results and point to an output folder

# Create the results directory
sinkdir <- 'Results'
dir.create(path = 'Results', showWarnings = F)

# Load the production set of results
output <- extract_output('Output')
output <- subset(output, grepl('_p$',Statistic))
output$Method <- gsub(x=output$Statistic, pattern='_p', repl='')

###################################################################################### ERROR RATES
# calculate the error rates
ER <- dcast(output, Method + Scenario + pSVS + decline ~ ., fun=error_rate, value.var='value')
names(ER)[ncol(ER)] <- 'ER' # default is to call it 'NA'

# get the sample sizes, just for a sense-check
ER$nreps <- as.numeric(acast(output, Method + Scenario + pSVS + decline ~ ., fun=length, value.var='value'))

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
print(ggV %+% subset(ER, decline==0 & grepl('A', Scenario))) # All scenario A
dev.off()

png(paste0(sinkdir,'/validity1 ',datecode,'.png'), wi=1024, hei=500)
print(ggV) # All validity 
dev.off()

png(paste0(sinkdir,'/validity2 ',datecode,'.png'), wi=1024, hei=500)
print(ggV %+% subset(ER, decline==0 &  Scenario != 'C1_pShortLEven'
               & !Method %in% c('nRecords', 'Maes', 'VRhov'))) # all except the really bad methods
dev.off()

png(paste0(sinkdir,'/validityMM ',datecode,'.png'), wi=1024, hei=500)
print(ggV %+% subset(ER, decline==0 & grepl('MM', Method))) # just MM
dev.off()

png(paste0(sinkdir,'/validity3 ',datecode,'.png'), wi=1024, hei=500)
print(ggV %+% subset(ER, decline==0 & Scenario != 'C1_pShortLEven'
               & !Method %in% c('nRecords', 'Maes')
               & !grepl('VR', Method)
               & !grepl('LL', Method) 
               & !grepl('A', Scenario)
               & !grepl('MMber', Method)
               )) # the good methods
dev.off()


# What happens just under B3
png(paste0(sinkdir,'/validityB3 ',datecode,'.png'), wi=500, hei=500)
qp <- qplot(data=subset(ER, decline==0 & grepl('B3f', Scenario) & Method %in% c('VRsimple', 'Occupancy', 'MMbin0', 'LL')),
      x=pSVS, y=ER, geom='line', col=Method, size=I(1))         
print(qp + geom_hline(yintercept =0.05))
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
print(ggP %+% subset(ER, grepl('A', Scenario) 
                     & !is.na(T2))) # All scenario A
dev.off()

png(paste0(sinkdir,'/Power1 ',datecode,'.png'), wi=1024, hei=500)
print(ggP) # All  
dev.off()

png(paste0(sinkdir,'/Power2 ',datecode,'.png'), wi=1024, hei=500)
print(ggP %+% subset(ER, Scenario != 'C1_pShortLEven'
                     & !Method %in% c('nRecords', 'Maes', 'VRhov')
                     & !is.na(T2))) # all except the really bad methods
dev.off()

png(paste0(sinkdir,'/PowerMM ',datecode,'.png'), wi=1024, hei=500)
print(ggP %+% subset(ER, grepl('MM', Method)
                     & !is.na(T2))) # just MM
dev.off()

png(paste0(sinkdir,'/Power3 ',datecode,'.png'), wi=1024, hei=500)
print(ggP %+% subset(ER, Scenario != 'C1_pShortLEven'
                     & !Method %in% c('nRecords', 'Maes')
                     & !grepl('VR', Method)
                     & !grepl('LL', Method) 
                     & !grepl('A', Scenario)
                     & !grepl('MMber', Method)
                     & !is.na(T2)))# the good methods
dev.off()




