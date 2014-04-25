# Nick Isaac 7 May 2013
# Extract and explore the simulation results
Explore_results <- function(readPath = 'Output', writePath = 'Results', plot = FALSE){

library(ggplot2)
library(reshape2)

homedir <- getwd()

source('extract_collate_stats.R')

datecode <- format(Sys.Date(), '%y%m%d')
outdir <- writePath

###################### Load the results and point to an output folder

# Create the results directory
sinkdir <- 'Results'
dir.create(path = 'Results', showWarnings = F)

# Load the production set of results
output_master <- extract_output(readPath)
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
  c('nSites', 'Naive', 1, T),
  c('Telfer', 'Telfer', 2, T),
  c('Frescalo1tp','Frescalo_Y', 4, T),
  c('Frescalo2tp','Frescalo_P', 3, T),
  c('VRsimple', 'ReportingRate', 5, T),
  c('LLsimple', 'RR+LL', 7, T),
  c('RR_SS', 'RR+SF', 6, T),
  c('MMbin0', 'RR+Site', 9, T),
  # two component models:
  c('LLmm', 'RR+LL+Site', 10, F),# in original ms. maybe include?
  c('mmSS', 'RR+SF+Site', 11, F),
  c('LLSS', 'RR+SF+LL', 12, F),
  
  # three component model
  c('LLmmSS', 'RR+SF+LL+Site', 18, T), # similar to WSS, but uses LL instead of visit threshold
  
  #Occupancy Models
  c('Occ+Simple','OccDetSimple', 20, T),
  c('Occ+LL+Site','Occ+LL+Site', 21, F), # virtually identical to the old 'Occupancy+Site' option
  c('Occ+SS+Site', 'Occ+SF+Site', 22, F),
  c('Occ+SS+LL', 'OD+SF+LL', 24, F),
  
  # site selection with three years
  c('RR_SS3', 'RR+SF3', 8, F),
  c('LLmmSS3', 'RR+LL+SF3+Site', 19, F), # similar to WSS, but uses LL instead of visit threshold
  c('Occ+SS3+Site', 'Occ+SF3+Site', 23, F),
  c('Occ+SS3+LL', 'Occ+LL+SF3', 25, F),
  
  c('Occ+Full', 'OD+SF+LL+Site', 26, T),
  
  #c('Occupancy','Occupancy', 12), # superseded
  #c('Occupancy+Site','Occupancy+Site', 13), # superseded
  c('MMbin2sp', 'WSS_2', 31, F), # to be dropped from resubmission & replaced with "RR+LL+Sel+Site"
  c('MMbin4sp', 'WSS_4', 32, F), # to be dropped from resubmission & replaced with "RR+LL+Sel+Site"
  c('Maes', 'RDC', 99, F)
  # plr, LLfs, Maes5, nSite, VRhov, nRecords: All in the output but ignored here
))
names(Methods_convert) <- c('sims', 'ms', 'id', 'include')
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

dim(ER) # 308   5

# calculate ER as an prdered factor
ER$intensity <- factor(ER$pSVS, levels=c(0.05,0.07,0.1), labels=c('low','med','high'), ordered=T)


##################################################################################### TYPE II ERROR RATES
# calculate the sample sizes & error rates
ER2 <- dcast(subset(output, sig_right_dir), Method + Scenario + pSVS ~ ., fun=length, value.var='p_val')
names(ER2)[ncol(ER2)] <- 'Nsig_correct' # default is to call it 'NA'

nrow(ER2) # <546, i.e. some combos have NO POWER

# merge
ER <- merge(ER, ER2, all=T)
subset(ER, is.na(Nsig_correct)) # no obvious pattern


################################# sample sizes

# get the sample sizes, just for a sense-check
ER$nrepsV <- as.numeric(acast(subset(output, decline==0), Method + Scenario + pSVS ~ ., fun=length, value.var='p_val'))

#ER$nrepsP <- as.numeric(acast(subset(output, decline==0.3), Method + Scenario + pSVS ~ ., fun=length, value.var='p_val'))
# some combos have no power so the above doesn't work
nrepsP <- dcast(subset(output, decline==0.3), Method + Scenario + pSVS ~ ., fun=length, value.var='p_val')
ER <- merge(ER, nrepsP, all.x=T)
names(ER)[ncol(ER)] <- 'nrepsP'

dcast(subset(ER, intensity=='med'), Method~Scenario, value.var='nrepsV')


####################################### type 2 error rate
ER$T2 <- 1 - with(ER, Nsig_correct/nrepsP)

# some of the rows have no power (i.e. no significant correlations in the right direction)
i <- is.na(ER$Nsig_correct)
ER$T2[i] <- 1

nrow(ER)  # 546 = 33 methods, 7 scenarios, 3 levels of intensity

# summarise T2 of all methods
acast(ER, Method ~ ., mean, value.var='T2')


###################################################################################### DEFINE POWER

# We define power here
ER$Power <- with(ER, 1 - T2 - T1)
ER$Power[ER$Power < 0] <- 0

######################################################################## MEDIAN TRENDS

# med_trends_v <- acast(
#   subset(output, decline==0 & pSVS==0.05 & Method %in% Methods_convert$ms & Scenario %in% with(Scenario_convert,ms[id<90])), 
#   Method ~ Scenario, value.var='trend', median, na.rm=T)
# 
# signif(med_trends_v, 3)

#good_methods <- c('WSS_2', 'WSS_4', 'RDC', 'Frescalo_Y', 'Frescalo_P', 'Telfer', 'Occupancy', 'Occupancy+Site')
#good_methods <- subset(Methods_convert, id < 30)$ms
#good_methods <- subset(Methods_convert, id < 19)$ms
#good_methods <- subset(Methods_convert, grepl('SS', ms))$ms
#good_methods <- subset(Methods_convert, grepl('Occ', ms))$ms
good_methods <- subset(Methods_convert, include==T)$ms

######################################################################## SHORTCUT
#load(file=paste0(outdir, '/ER_',datecode,'.rData')) ## ->ER
######################################################################## PREPARE TO PLOT

# Setup x axis labes
ticks <- data.frame(t = c(0.05,0.07,0.1),
                    l = c('L', 'M', 'H'))

# Setup colours
if(length(good_methods) == 13){
  colours <- c("#F8766D", "#CD9600", 
               "#579438", "#82DB56", # greens
               "#43A39A", "#6AF7E9", "#5E95C4", "#A6D1F7", # blues
               "#8258CC", "#9C6E9B", "#E8A5E7", "#AD3687", "#FF61CC") # purple-pinks
} else if(length(good_methods) == 11) {
  colours <- c("#F8766D", "#CD9600", "#579438", "#82DB56",
               "#43A39A", "#6AF7E9", "#5E95C4", "#A6D1F7",
               "#8258CC", # purple on screen, blue on printed page
               "#AD3687", "#FF61CC")
} else {
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  colours <- gg_color_hue(length(good_methods))
}
names(colours) <- sort(good_methods)
colScale <- scale_colour_manual(name = "Method",values = colours)


######################################################################## BARPLOT
# build plot
ggV <- ggplot(data=ER[ER$intensity == 'med',], aes(x=Method, y = T1))
ggV <- ggV + facet_wrap(~Scenario, nrow=2)
ggV <- ggV + geom_bar(aes(fill = Method),stat="identity")
ggV <- ggV + scale_fill_manual(values = colours)
#ggV <- ggV + scale_y_log10(limits = c(1,1000), expand=c(0,0), breaks = c(1, 10, 100, 1000),labels = c('0.001', '0.01','0.1','1'))
ggV <- ggV + scale_y_sqrt(breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1))
#ggV <- ggV + scale_y_sqrt()
ggV <- ggV + geom_hline(yintercept = 0.05, lwd = 0.5,lty = 1)
ggV <- ggV + geom_hline(yintercept = 0.1, lwd = 0.5,lty = 2)
ggV <- ggV + ylab('Type I error rate') + xlab('Method')
ggV <- ggV + theme_bw() 
ggV <- ggV + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank())

#ggV # all methods & scenarios
png(paste0(outdir,'/Validity ',datecode,'.png'), wi=1024, hei=500)
print(ggV %+% subset(ER[ER$intensity == 'med',], Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))
dev.off()
if (plot) print(ggV %+% subset(ER[ER$intensity == 'med',], Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))
######################################################################## BAR PLOTS POWER

ER$Power[ER$Power == 0.000] <- 0.005

# build plot
ggP <- ggplot(data=ER[ER$intensity == 'med',], aes(x=Method, y = Power))
ggP <- ggP + facet_wrap(~Scenario, nrow=2)
ggP <- ggP + geom_bar(aes(fill = Method),stat="identity")
ggP <- ggP + scale_fill_manual(values = colours)
#ggP <- ggP + ylab('Power') #+ xlab('Recording Intensity')
ggP <- ggP + theme_bw() 
ggP <- ggP + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank())

#ggV # all methods & scenarios
png(paste0(outdir,'/Power ',datecode,'.png'), wi=1024, hei=500)
print(ggP %+% subset(ER[ER$intensity == 'med',], Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))
dev.off()
if(plot) print(ggP %+% subset(ER[ER$intensity == 'med',], Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))

######################################################## trend estimates

#mean_trends <- acast(subset(trends), Method~Scenario~pSVS, mean)
#round(mean_trends[,,],3)


######################################################################## OUTPUT TABLES

for(level in unique(ER$pSVS)){ # pSVS at which to report the results
  
  # simple output table of medium recording intensity
  (Vtable <- dcast(subset(ER, pSVS==level & Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])), 
                   Method ~ Scenario, value.var='T1'))
  write.csv(Vtable, file=paste0(outdir,'/table_V_',level,'.csv'))
  
  # simple power table of medium recording intensity
  (Ptable <- dcast(subset(ER, pSVS==level & Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])), 
                   Method ~ Scenario, value.var='Power'))
  write.csv(Ptable, file=paste0(outdir,'/table_P_',level,'.csv'))
  
}

######################################################## SAVE

#write.csv(ER, file=paste0(outdir, '/ER_',datecode,'.csv'))
save(ER, file=paste0(outdir, '/ER_',datecode,'.rData'))

######################################################################## SI Figures

ggV <- ggplot(data=ER) 
ggV <- ggV + scale_x_continuous(breaks=c(ticks$t), 
                                labels=as.character(ticks$l))
ggV <- ggV + geom_line(aes(x=pSVS, y=T1, col=Method), size=1)
ggV <- ggV + scale_y_sqrt(breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1))
ggV <- ggV + geom_hline(yintercept = 0.05, lwd = 0.5,lty = 1)
ggV <- ggV + geom_hline(yintercept = 0.1, lwd = 0.5,lty = 2)
ggV <- ggV + ylab('Type I error rate') + xlab('Recording Intensity')
#ggV <- ggV + geom_point(aes(col=Method, size=3))
ggV <- ggV + colScale + theme_bw() 
ggV <- ggV + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggV <- ggV + facet_wrap(~Scenario, nrow=1) 

#ggV # all methods & scenarios
ggV %+% subset(ER, Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90]))

# now write to a file
png(paste0(outdir,'/validity_SI.png'), wi=1024, hei=500)
print(ggV %+% subset(ER, Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))
dev.off()
if(plot) print(ggV %+% subset(ER, Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))


ggP <- ggplot(data=ER, aes(x=pSVS, y=Power, gr=Method)) # try replacing with intensity
ggP <- ggP + scale_x_continuous(breaks=c(ticks$t), 
                                labels=as.character(ticks$l))
ggP <- ggP + geom_line(aes(col=Method), size=1)
#ggP <- ggP + scale_y_log10() + geom_hline(yintercept =0.05)
ggP <- ggP + ylab('Power') + xlab('Recording Intensity')
#ggP <- ggP + geom_point(aes(col=Method, size=3))
#ggP # all methods & scenarios#
ggP <- ggP + colScale
ggP <- ggP + theme_bw()
ggP <- ggP + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggP <- ggP + facet_wrap(~Scenario, nrow=1)
ggP %+% subset(ER, Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90]))

# now write to a file
png(paste0(outdir,'/power_SI.png'), wi=1024, hei=500)
print(ggP %+% subset(ER, Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))
dev.off()
if(plot) print(ggP %+% subset(ER, Method %in% good_methods & Scenario %in% with(Scenario_convert,ms[id<90])))

######################################################################## Figures 2

# now write to a file
png(paste0(outdir,'/figure 3.png'), wi=500, hei=500)
print(ggP %+% subset(ER, Method %in% good_methods & Scenario == 'Control') + ylim(c(0,1)))
dev.off()
if(plot) print(ggP %+% subset(ER, Method %in% good_methods & Scenario == 'Control') + ylim(c(0,1)))
}