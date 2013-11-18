model {                       

# Priors 

# state model priors

for(t in 1:nyear){
    a[t] ~ dunif(-10,10)   
}                 

for (i in 1:nsite) {
   eta[i] ~ dnorm(mu2, tau2)       # extra random site-effect on occupancy
} 

mu2 ~ dnorm(0, 0.001)
tau2 <- 1/(sigma2 * sigma2)
sigma2 ~ dunif(0, 5)


# observation model priors 
for (t in 1:nyear) {
    alpha.p[t] ~ dnorm(mu.lp, tau.lp)            # p random year
}

mu.lp ~ dnorm(0, 0.01)                         
tau.lp <- 1 / (sd.lp * sd.lp)                 
sd.lp ~ dunif(0, 5)   
                           
dtype2.p ~ dunif(dtype2p_min,dtype2p_max) 

# State model

for (i in 1:nsite){ 
     for (t in 1:nyear){   
      z[i,t] ~ dbern(muZ[i,t]) # True occupancy z at site i
      logit(muZ[i,t])<- a[t] + eta[i] # plus random site effect
   } 
}   


# Observation model 

for (t in 1:nyear){
     for(i in 1:nsite){
         for(j in 1:nvisit) {
            Py[i,j,t]<- z[i,t]*p[i,j,t]
  logit(p[i,j,t]) <- alpha.p[t] + dtype2.p*DATATYPE2[i,j,t] 
            M[i,j,t] ~ dbern(Py[i,j,t])  

Presi[i,j,t] <- abs(M[i,j,t]-p[i,j,t])
y.new[i,j,t] ~ dbern(Py[i,j,t]) 
Presi.new[i,j,t] <- abs(y.new[i,j,t]-p[i,j,t])

         } 
     }
 } 

# Bayesian Goodness-of-Fit
fit<-sum(Presi[,,])
fit.new<- sum(Presi.new[,,])

# Derived parameters state model

# Finite sample occupancy
for (t in 1:nyear) {  
    psi.fs[t] <- sum(z[1:nsite,t])/nsite
} 

# Overall trend in occpuancy
sumY <- sum(psi.fs[1:nyear])
for (t in 1:nyear) {
   sumxy[t] <- psi.fs[t]*t
}
sumXY <- sum(sumxy[1:nyear])
regres.psi <- (sumXY - ((sumX*sumY)/nyear))/(sumX2 - ((sumX*sumX)/nyear))

# Derived parameters observation model
for (t in 1:nyear) {          
   pdet.alpha[t] <- exp(alpha.p[t])/(1 + exp(alpha.p[t])) 
   pdet.d2[t] <- exp(alpha.p[t]+dtype2.p)/(1 + exp(alpha.p[t]+dtype2.p))
}

# overall trend in pdet.alpha
sumYpdet <- sum(pdet.alpha[1:nyear])  
for (t in 1:nyear) {          
      sumxypdet[t] <- pdet.alpha[t]*t
}
sumXYpdet <- sum(sumxypdet[1:nyear])
regres.pdet <- (sumXYpdet - ((sumX*sumYpdet)/nyear))/(sumX2 - ((sumX*sumX)/nyear))

# end of model formulation
}  				