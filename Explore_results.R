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
output_master <- extract_output('Output')
output_master$Statistic <- gsub(output_master$Statistic, pattern='Telfer$', repl='Telfer_trend')
output_master$Statistic <- gsub(output_master$Statistic, pattern='plr$', repl='plr_trend')
output <- subset(output_master, grepl('_p$',Statistic))
output$Method <- gsub(x=output$Statistic, pattern='_p', repl='')

trends <- subset(output_master, grepl('_trend$',output_master$Statistic))
trends$Method <- gsub(x=trends$Statistic, pattern='_trend', repl='')

# stick them together
merge_cols <- c('Method', 'Scenario',  'repnum', 'decline', 'pSVS')
output <- merge(output[,-1], trends[,c(merge_cols, 'value')], by=merge_cols)

# rename the key data columns
names(output) <- gsub(names(output), patt='value.x', repl='p_val')
names(output) <- gsub(names(output), patt='value.y', repl='trend')

# is the trend signficant & in the right direction?
output$sig_right_dir <- with(output, p_val <= 0.05 & trend * decline < 0) # returns FALSE for all validity tests

# We're using different names in the Paper than the code
Methods_convert <- as.data.frame(rbind(
  c('LL', 'ListLength', 7),
  c('LLmm', 'ListLength+Site', 8),
  c('nSites', 'Naive', 1),
  c('VRsimple', 'ReportingRate', 5),
  c('MMbin0', 'Reporting+Site', 6),
  c('MMbin2sp', 'WSS_2', 3),
  c('MMbin4sp', 'WSS_4', 4),
  c('Telfer', 'Telfer', 8),
  c('Frescalo1tp','Frescalo_Y', 11),
  c('Frescalo2tp','Frescalo_P', 10),
  c('Occupancy','Occupancy', 12),
  c('Occupancy+Site','Occupancy+Site', 13),
  c('Maes', 'RDC', 2)
))
names(Methods_convert) <- c('sims', 'ms', 'id')
Methods_convert$id <- as.numeric(as.character(Methods_convert$id))

# which methods are in the sims but excluded from above
unused_methods <- setdiff(output$Method, Methods_convert$sims)

# reorder the factor levels to the order above
Methods_convert$ms <- with(Methods_convert, factor(ms, ordered=T, levels=ms[order(id)]))
# now append the unused methods as levels
levels(Methods_convert$ms) <- c(levels(Methods_convert$ms), unused_methods) # works a treat

output$Method_sims <- output$Method
output$Method <- Methods_convert$ms[match(output$Method, Methods_convert$sims)]

i <- is.na(output$Method)
output$Method[i] <- output$Method_sims[i]

###################################################################################### RENAME SCENARIOS

Scenario_convert <- as.data.frame(rbind(
  c('A_EvenRcrdng', 'Control', 1),
  c('B2_IncrnVisit', 'MoreVisits', 2),
  c('B3f_IncrnVBiasFc', 'MoreVisits+Bias', 3),
  c('C1_pShortLEven', 'C1_pShortLEven', 99),
  c('C2_pShortLIncr', 'LessEffortPerVisit', 4),
  c('D2_SelectvIncr', 'MoreDetectable', 5),
  c('F_NfDecline', 'NonFocalDeclines', 9)
))
names(Scenario_convert) <- c('sims', 'ms', 'id')
Scenario_convert$id <- as.numeric(as.character(Scenario_convert$id))

# reorder the factor levels
Scenario_convert$ms <- with(Scenario_convert, factor(ms, ordered=T, levels=ms[order(id)]))

output$Scenario <- Scenario_convert$ms[match(output$Scenario, Scenario_convert$sims)]

###################################################################################### TYPE I ERROR RATES
# calculate the error rates
ER <- dcast(subset(output, decline==0), Method + Scenario + pSVS ~ ., fun=error_rate, value.var='p_val')
names(ER)[ncol(ER)] <- 'T1' # default is to call it 'NA'

# NB under some scenarios, WSS_4 has p_val=NA. This is where there's not enough data to fit the model

# get the sample sizes, just for a sense-check
ER$nreps <- as.numeric(acast(subset(output, decline==0), Method + Scenario + pSVS ~ ., fun=length, value.var='p_val'))

dim(ER) # 378   5


##################################################################################### TYPE II ERROR RATES
# calculate the error rates
ER2 <- dcast(subset(output, sig_right_dir), Method + Scenario + pSVS ~ ., fun=length, value.var='p_val')
names(ER2)[ncol(ER2)] <- 'Nsig_correct' # default is to call it 'NA'

nrow(ER2) # 365 i.e. 13 combos have NO POWER

# merge
ER <- merge(ER, ER2, all=T)
subset(ER, is.na(Nsig_correct)) # no obvious pattern

ER$T2 <- 1 - with(ER, Nsig_correct/nreps)

# some of the rows have no power (i.e. no significant correlations in the right direction)
i <- is.na(ER$Nsig_correct)
ER$T2[i] <- 1

# calculate ER as an prdered factor
ER$intensity <- factor(ER$pSVS, levels=c(0.05,0.07,0.1), labels=c('low','med','high'), ordered=T)

nrow(ER)  # 378 = 18 methods, 7 scenarios, 3 levels of intensity

# summarise T2 of all methods
acast(ER, Method ~ ., mean, value.var='T2')

###################################################################################### DEFINE POWER

# We define power here
ER$Power <- with(ER, 1 - T2 - T1)
ER$Power[ER$Power < 0] <- 0

######################################################################## OUTPUT TABLES

level <- 0.07 # pSVS at which to report the results

# simple output table of medium recording intensity
(table3 <- dcast(subset(ER, pSVS==level & Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90])), 
                 Method ~ Scenario, value.var='T1'))
write.csv(table3, file='Results/Table_S1.csv')

# simple power table of medium recording intensity
(table4 <- dcast(subset(ER, pSVS==level & Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90])), 
                 Method ~ Scenario, value.var='Power'))
write.csv(table4, file='Results/Table_S2.csv')


######################################################################## MEDIAN TRENDS

med_trends_v <- acast(
  subset(output, decline==0 & pSVS==0.05 & Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90])), 
  Method ~ Scenario, value.var='trend', median, na.rm=T)

signif(med_trends_v, 3)

good_methods <- c('WSS_2', 'WSS_4', 'RDC', 'Frescalo_Y', 'Frescalo_P', 'Telfer', 'Occupancy', 'Occupancy+Site')

######################################################################## PREPARE TO PLOT

#code by Tom August
# Setup x axis labes
ticks <- data.frame(t = c(0.05,0.07,0.1),
                    l = c('L', 'M', 'H'))

# Set up colours
colours <- c("#F8766D", "#CD9600", "#579438", "#82DB56", "#43A39A",
             "#6AF7E9", "#5E95C4", "#A6D1F7", "#8258CC", "#9C6E9B",
             "#E8A5E7", "#AD3687", "#FF61CC")
names(colours) <- sort(Methods_convert$ms)
colScale <- scale_colour_manual(name = "Method",values = colours)

########################################################### PLOTTING VALIDITY WITH INTENSITY
ggV <- ggplot(data=ER) 
ggV <- ggV + scale_x_continuous(breaks=c(ticks$t), 
                                labels=as.character(ticks$l))
ggV <- ggV + geom_line(aes(x=pSVS, y=T1, col=Method), size=1)
ggV <- ggV + scale_y_sqrt(breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1))
ggV <- ggV + geom_hline(yintercept = 0.05, lwd = 0.5,lty = 1)
ggV <- ggV + geom_hline(yintercept = 0.1, lwd = 0.5,lty = 2)
ggV <- ggV + ylab('Type I error rate') + xlab('Recording Intensity')
ggV <- ggV + colScale + theme_bw() 
ggV <- ggV + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggV <- ggV + facet_wrap(~Scenario, nrow=1) 


# MAKE SOME PLOTS OF VALIDITY
png(filename='Results/Figure_S1.png', wi=1024, hei=500)
print(ggV %+% subset(ER, Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90])))
dev.off()

########################################################### PLOTTING VALIDITY WITHOUT INTENSITY
# build plot
ggV <- ggplot(data=ER[ER$intensity == 'med',], aes(x=Method, y = T1))
ggV <- ggV + facet_wrap(~Scenario, nrow=2)
ggV <- ggV + geom_bar(aes(fill = Method),stat="identity")
ggV <- ggV + scale_fill_manual(values = colours)
ggV <- ggV + scale_y_sqrt(breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1))
ggV <- ggV + geom_hline(yintercept = 0.05, lwd = 0.5,lty = 1)
ggV <- ggV + geom_hline(yintercept = 0.1, lwd = 0.5,lty = 2)
ggV <- ggV + ylab('Type I error rate') + xlab('Method')
ggV <- ggV + theme_bw() 
ggV <- ggV + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank())

# MAKE SOME PLOTS OF VALIDITY
png(filename='Results/Figure_1.png', wi=1024, hei=500)
print(ggV %+% subset(ER[ER$intensity == 'med',], Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90])))
dev.off()


################################################################## NOW REPEAT FOR POWER

# build the ggplot object
ggP <- ggplot(data=ER, aes(x=pSVS, y=Power, gr=Method)) 
ggP <- ggP + scale_x_continuous(breaks=c(ticks$t), 
                                labels=as.character(ticks$l))
ggP <- ggP + geom_line(aes(col=Method), size=1)
ggP <- ggP + ylab('Power') + xlab('Recording Intensity')
ggP <- ggP + colScale
ggP <- ggP + theme_bw()
ggP <- ggP + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggP <- ggP + facet_wrap(~Scenario, nrow=1)
ggP %+% subset(ER, Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90]))


# MAKE SOME PLOTS OF POWER

png(filename='Results/Figure_2.png', wi=500, hei=250)
print(ggP %+% subset(ER, Scenario=='Control' & Method %in% good_methods))
dev.off()

png(filename='Results/Figure_S2.png', wi=1024, hei=500)
print(ggP %+% subset(ER, Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90])))
dev.off()

######################################################################## BAR PLOTS POWER

ER$Power[ER$Power == 0.000] <- 0.005

# build plot
ggP <- ggplot(data=ER[ER$intensity == 'med',], aes(x=Method, y = Power))
ggP <- ggP + facet_wrap(~Scenario, nrow=2)
ggP <- ggP + geom_bar(aes(fill = Method),stat="identity")
ggP <- ggP + scale_fill_manual(values = colours)
ggP <- ggP + theme_bw() 
ggP <- ggP + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank())

#ggV # all methods & scenarios
ggP %+% subset(ER[ER$intensity == 'med',], Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90]))

# MAKE SOME PLOTS OF POWER
png(filename='Results/Figure_3.png', wi=1024, hei=500)
print(ggP %+% subset(ER[ER$intensity == 'med',], Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))
dev.off()