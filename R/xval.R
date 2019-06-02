xvalPoly <- function(reps=5, 
                     mean_loss=function(y, yhat){mean(abs(y - yhat))},
                     stat=median, pTesting=0.1, seed=NULL, ...){
  
  args <- list(...)
  if(!('xy' %in% names(args)))
    stop('Please provide xy and deg and any other argument required by polyFit().')
  if(!('deg' %in% names(args)))
    stop('Please provide xy and deg and any other argument required by polyFit().')
  args$noisy <- if(is.null(args$noisy)) TRUE else args$noisy
    
  n <- nrow(args$xy)
  
  if(!is.null(seed)) set.seed(seed)
    
  mean_test_loss <- c() # out of sample average loss
  
  if(args$noisy) cat("##############################\nBeginning cross-validations for deg", 
                     args$deg, "and maxInteractDeg", 
                     if('maxInteractDeg' %in% names(args)) args$maxInteractDeg else args$deg,
                     "\n\n")
  
  for(i in 1:reps){
    
    test <- sample(c(TRUE, FALSE), n, replace=TRUE, prob = c(pTesting, 1 - pTesting))
    ytest <- args$xy[test, ncol(args$xy)]
    
    trained <- polyFit(xy = args$xy[!test, ], 
                       deg = args$deg, 
                       maxInteractDeg = if('maxInteractDeg' %in% names(args)) args$maxInteractDeg else args$deg, 
                       use = if('use' %in% names(args)) args$use else 'lm', 
                       pcaMethod = if('pcaMethod' %in% names(args)) args$pcaMethod else NULL,
                       pcaLocation = if('pcaLocation' %in% names(args)) args$pcaLocation else 'front', 
                       pcaPortion = if('pcaPortion' %in% names(args)) args$pcaPortion else 0.9, 
                       glmMethod = if('glmMethod' %in% names(args)) args$glmMethod else'one', 
                       return_xy = if('return_xy' %in% names(args)) args$return_xy else FALSE, 
                       returnPoly = if('returnPoly' %in% names(args)) args$returnPoly else FALSE, 
                       noisy = args$noisy
                       )
    
    yhat <- predict(trained, newdata=args$xy[test, ])
    mean_test_loss[i] <- mean_loss(ytest, yhat)
    if(args$noisy)
      cat("\n\n\nFinished iteration", i, "\nMean loss on test data:", mean_test_loss[i], "\n")
    
  }
  
  typical_loss <- stat(mean_test_loss)
  if(args$noisy)
    cat("\n\n\n\nFinished.\n\nTypical loss (as defined by user inputed stat) over repetitions:", typical_loss, "\n\n")
  return(typical_loss)
  
    
}