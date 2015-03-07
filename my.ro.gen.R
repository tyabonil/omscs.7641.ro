library(nnet)
library(GenSA)

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

nnet.rhc <- function(rhc.d=0.01, rhc.min=-0.5, rhc.max=0.5, rhc.n=3, rhc.test, rhc.testIdx, rhc.seed=1, ...) {
  dots <- list(...)
  rhc.scores <- NULL
  j <- 1
  rhc.best.score <- 0
  
  rhc.keepgoing  <- TRUE  
  while(all(rhc.keepgoing)==TRUE) {
    rhc.keepgoing <- FALSE
    set.seed(rhc.seed*j)
    #print(c(j, rhc.best.score))
    j <- j + 1 
    rhc.initial <- runif(n=rhc.n, min=rhc.min, max=rhc.max)
    rhc.neighbors <- rhc.nb(val=rhc.initial, d=rhc.d)
    i <- 1
    rhc.best.value <- 0
    rhc.best.fit <- NULL
    rhc.best.results <- NULL
    for(i in 1:nrow(rhc.neighbors)) {
    #for(i in 1:1) {
      rhc.test.value <- rhc.neighbors[i, ]
      rhc.test.fit <- nnet(..., Wts=rhc.test.value)
      rhc.test.results <- predict(rhc.test.fit, newdata=rhc.test, type="class")
      rhc.test.score <- sum((rhc.test.results == rhc.testIdx) == TRUE)      
      rhc.scores <- rbind(rhc.scores, rhc.test.score)
      if(rhc.test.score > rhc.best.score) {
        rhc.keepgoing <- TRUE
        rhc.best.value <- rhc.test.value
        rhc.best.fit <- rhc.test.fit
        rhc.best.results <- rhc.test.results
        rhc.best.score <- rhc.test.score
        
      } 
    }
  }
  list(fit=rhc.best.fit, results=rhc.best.results, score=rhc.best.score, scores=rhc.scores) 
}

## sa.temp is initial temperature
## sa.test is the test set
## sa.testIdx is the test results expected
## sa.seed is the start seed
## sa.d is the absolute value of the upper and lower bounds to the weights
## sa.n is the number of weights

nnet.sa <- function(sa.temp, sa.test, sa.testIdx, sa.seed=1, sa.d=100, sa.n=19,...) {
  set.seed(sa.seed)
  lower <- as.numeric(-1.*rep(sa.d, sa.n))
  upper <- as.numeric(-1.*lower)
  fit <- function(fit.wts, fit.test=sa.test, fit.testIdx=sa.testIdx) {
    temp.fit <- nnet(..., Wts=fit.wts)
    temp.predict <- predict(temp.fit, newdata=fit.test, type="class")
    temp.score <- sum((temp.predict == fit.testIdx) == FALSE)
    as.numeric(temp.score)
  }
  GenSA(par=NULL, fn=fit, lower=lower, upper=upper,
        control=list(max.time=10, temperature=nrow(sa.test)))
}

