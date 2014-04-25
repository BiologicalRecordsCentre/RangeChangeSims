# Runs a JAGS test as a way to check that it is installed
jagsTest<-function(){
  
  x<-try({
    model.file <- system.file(package="R2jags", "model", "schools.txt")
    J <- 8.0
    y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
    sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6)
    jags.data <- list("y","sd","J")
    jags.params <- c("mu","sigma","theta")
    jags.inits <- function(){
      list("mu"=rnorm(1),"sigma"=runif(1),"theta"=rnorm(J))
    }
    jagsfit <- jags(data=list("y","sd","J"), inits=jags.inits, jags.params,
                    n.iter=10, model.file=model.file)
  }, silent = TRUE)
  
  if(class(x) == 'try-error'){
    stop('The test JAGS script did not run. Please ensure you have installed JAGS')
  } else if(class(x)){
    cat('JAGS test successful/n')
  } else{
    warning(paste('JAGS test returned object of class', class(x), sep = ' '))
  }
}