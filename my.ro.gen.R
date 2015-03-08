library(nnet)
library(GenSA)
library(GA)

## This function creates a list of neighbours to test
## Distance function adds a value to each dimension

rhc.nb <- function(val=0, d=0.01) {
  nb.out <- val
  n <- length(val)
  max <- d # runif(1, 0, d)
  min <- -1*max
  #nb.out <- rbind(nb.out, nb.out + min)
  #nb.out <- rbind(nb.out, nb.out + max) 
  for(i in 1:n) {
    nb.left <- val
    nb.left[i] <- nb.left[i] - d
    nb.right <- val
    nb.right[i] <- nb.right[i] + d    
    nb.out <- rbind(nb.out, nb.left)
    nb.out <- rbind(nb.out, nb.right)    
  }
  nb.out
}
## rhc.test is the test set
## rhc.testIdx is the test results expected
## rhc.seed is the start seed
## rhc.min/max are the min and max values of the weights
## rhc.d is the distance threshold to evaluate a neighbour
## rhc.n is the number of weights
## rhc.iter is the number of iterations that we loop over


nnet.rhc <- function(rhc.d=0.01, rhc.min=-0.5, rhc.max=0.5, rhc.n=3, rhc.test, rhc.testIdx, rhc.seed=1, rhc.iter=100, ...) {
  t1 <- proc.time()
  j <- 1
  rhc.best.score <- 0
  
  rhc.keepgoing  <- TRUE
  for(k in 1:rhc.iter) {
    while(all(rhc.keepgoing)==TRUE) {
      rhc.keepgoing <- FALSE
      set.seed(rhc.seed*k)
      #print(c(j, rhc.best.score))
      j <- j + 1 
      rhc.initial <- runif(n=rhc.n, min=rhc.min, max=rhc.max)
      rhc.neighbors <- rhc.nb(val=rhc.initial, d=rhc.d)
      i <- 1
      rhc.best.value <- 0
      rhc.best.results <- NULL
      for(i in 1:nrow(rhc.neighbors)) {
        #for(i in 1:1) {
        rhc.test.value <- rhc.neighbors[i, ]
        rhc.test.fit <- nnet(..., Wts=rhc.test.value)
        rhc.test.results <- predict(rhc.test.fit, newdata=rhc.test, type="class")
        rhc.test.score <- sum((rhc.test.results == rhc.testIdx) == TRUE)      
        if(rhc.test.score > rhc.best.score) {
          rhc.keepgoing <- TRUE
          rhc.best.value <- rhc.test.value
          rhc.best.fit <- rhc.test.fit
          rhc.best.results <- rhc.test.results
          rhc.best.score <- rhc.test.score     
          rhc.best.seed <- rhc.seed*k
        } 
      }
    }  
  }
  
  rhc.best.score <- rhc.best.score/nrow(rhc.test)
  
  rhc.run.time <- proc.time() - t1
  list(best.model=rhc.best.fit, best.score=rhc.best.score, best.seed=rhc.best.seed, run.time=rhc.run.time) 
}

## sa.temp is initial temperature
## sa.test is the test set
## sa.testIdx is the test results expected
## sa.seed is the start seed
## sa.d is the absolute value of the upper and lower bounds to the weights
## sa.n is the number of weights
## sa.iter is the number of iterations that is passed to the SA function


nnet.sa <- function(sa.temp, sa.test, sa.testIdx, sa.seed=1, sa.d=100, sa.n=19, sa.iter=100, ...) {
  t1 <- proc.time()
  
  set.seed(sa.seed)
  lower <- as.numeric(-1.*rep(sa.d, sa.n))
  upper <- as.numeric(-1.*lower)
  fit <- function(fit.wts, fit.test=sa.test, fit.testIdx=sa.testIdx) {
    temp.fit <- nnet(..., Wts=fit.wts)
    temp.predict <- predict(temp.fit, newdata=fit.test, type="class")
    temp.score <- sum((temp.predict == fit.testIdx) == FALSE)
    as.numeric(temp.score)
  }
  best.fit <- GenSA(par=NULL, fn=fit, lower=lower, upper=upper,
        control=list(maxit=sa.iter, temperature=nrow(sa.test)))
  sa.best.fit <- nnet(..., Wts=best.fit$par)
  sa.best.score <- (nrow(sa.test) - best.fit$value)/nrow(sa.test)
  sa.best.seed <- sa.seed
  sa.run.time <- proc.time() - t1
  list(best.model=sa.best.fit, best.score=sa.best.score, best.seed=sa.best.seed, run.time=sa.run.time)  
}

## ga.test is the test set
## ga.testIdx is the test results expected
## ga.seed is the start seed
## ga.d is the absolute value of the upper and lower bounds to the weights
## ga.n is the number of weights
## ga.iter is the number of iterations that is passed to the GA function

nnet.ga <- function(ga.test, ga.testIdx, ga.seed=1, ga.d=100, ga.n=19, ga.iter=10, ...) {
  t1 <- proc.time()
  
  set.seed(ga.seed)
  lower <- as.numeric(-1.*rep(ga.d, ga.n))
  upper <- as.numeric(-1.*lower)
  fit <- function(fit.wts, fit.test=ga.test, fit.testIdx=ga.testIdx) {
    temp.fit <- nnet(..., Wts=fit.wts)
    temp.predict <- predict(temp.fit, newdata=fit.test, type="class")
    temp.score <- sum((temp.predict == fit.testIdx) == TRUE)
    as.numeric(temp.score)
  }
  best.fit <- ga(type="real-valued", fitness=fit, min=lower, max=upper, maxiter=ga.iter, monitor=NULL
        )
  ga.best.fit <- nnet(..., Wts=best.fit@solution[1,])
  ga.best.score <- max(best.fit@fitness)/nrow(ga.test)
  ga.best.seed <- ga.seed
  ga.run.time <- proc.time() - t1
  list(best.model=ga.best.fit, best.score=ga.best.score, best.seed=ga.best.seed, run.time=ga.run.time)  
}

