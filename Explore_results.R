# Nick Isaac 19-20 June 2013
# Extract and explore the simulation results
# We anwt to know if fixed versus variabel number of visits has an  
# this time we'll collate the output

###########

# Previously (130619) I compared the most recent results with 130502, which was a test run
# Actually I should compare them with the 'production runs' of 130509
# In that set, I *thought* that the number of visits was fixed, but actually it was variable

###########


library(ggplot2)
library(reshape2)
rm(list=objects())
homedir <- 'P:/NEC04273_SpeciesDistribution/Workfiles/Range change sims' # same location as 'Sim functions 13xxxx.r'
setwd(homedir)

source('extract_collate_stats.R')
source('Sim_functions')


sinkdir <- 'results/resuts 130619 stochastic=F'


nrow(output_stoch <- extract_output('results/results 130509')) #1155000
nrow(output_fixed <- extract_output('results/resuts 130619 stochastic=F')) #273000


# check the names match
all(names(output_stoch) == names(output_fixed)) # T

# I belived the 'stoch' set was created with a fixed number of visits per year. 
# In fact I made a coding error so it was actually variable (i.e. stochastic)
output_stoch$stochastic <- 1


# which methods have been used?
setdiff(output_stoch$Statistic, output_fixed$Statistic) # Frescalo & MM absent from the 'fixed' set

# restrict the output to just the columns of interest (i.e. the 'Methods')
nrow(output_fixed <- subset(output_fixed, grepl('_p$',Statistic))) #56000

# and do the same in the variable set
nrow(output_stoch <- subset(output_stoch, grepl('_p$',Statistic)
                            & Statistic %in% unique(output_fixed$Statistic)
                            #& pSVS %in% unique(output_fixed$pSVS)
                            )) #168000

# how big are they
with(output_stoch, table(decline, pSVS)) # 2 times bigger
with(output_fixed, table(decline, pSVS))

# how many reps per setup?
max(output_stoch$repnum) # 500!
max(output_fixed$repnum) # only 250!

     
output <- droplevels(rbind(output_fixed, output_stoch))
table(output$stochastic) #8000 F, 7000 T (mismatch is becuase later set includes Maes method)

output$Method <- gsub(x=output$Statistic, pattern='_p', repl='')



######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## 



# calculate the error rates
ER <- dcast(output, Method + Scenario + pSVS + stochastic + decline ~ ., fun=error_rate, value.var='value')
names(ER)[ncol(ER)] <- 'ER' # default is to call it 'NA'

# get the sample sizes, just for a sense-check
(ERN <- dcast(output, Method + Scenario + pSVS + stochastic + decline ~ ., fun=length, value.var='value'))

ER$combo <- with(ER, paste(stochastic,Method))

# build the ggplot object
ggV <- ggplot(data=subset(ER, decline==0), aes(x=pSVS, y=ER, gr=combo))
ggV <- ggV + geom_line(aes(col=Method, gr=combo, lty=factor(stochastic)), size=1)
ggV <- ggV + scale_y_log10() + geom_hline(yintercept =0.05)
ggV <- ggV + ylab('Type I error rate') + xlab('recording intensity')
ggV <- ggV + geom_point(aes(col=Method, shape=factor(stochastic)), size=3)
ggV <- ggV + facet_wrap(~Scenario, nrow=2)
ggV


png(paste0(sinkdir,'/ER ~stoch 130620.png'), wi=1024, hei=500)
ggV %+% subset(ER, decline==0 & Scenario != 'C1_pShortLEven'
               & !Method %in% c('nRecords', 'Maes', 'VRhov'))
dev.off()