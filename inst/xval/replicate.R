rep_xvalPoly <- function(reps = 5, stat = median, nHoldout = NULL, yCol = NULL, 
                         seed = NULL, estimateLowerDegrees = TRUE, ...){
  
  args <- list(...)
  if(!('xy' %in% names(args)))
    stop('Please provide xy and deg and any other argument for polyFit() as named arguments. Example:\n\nrep_xvalPoly(xy=iris, deg=2, use="glm", glmMethod="multlog")')
  if(!('deg' %in% names(args)))
    stop('Please provide xy and deg and any other argument for polyFit() as named arguments. Example:\n\nrep_xvalPoly(xy=iris, deg=2, use="glm", glmMethod="multlog")')
  args$noisy <- if(is.null(args$noisy)) TRUE else args$noisy
  
  xy <- args$xy
  
  n <- nrow(xy)
  nHoldout <- min(10000, round(0.2*n))
  
  if(!is.null(yCol)){
    xy <- cbind(xy[,-yCol], xy[,yCol])
  }else{
    yCol <- ncol(xy)
  }
    
  if(length(unique(xy[,yCol])) == 2 || is.character(xy[,yCol]))
    xy[,yCol] <- as.factor(xy[,yCol])
  
  if(!is.null(seed)) set.seed(seed)
  
  mean_test_loss <- matrix(nrow = reps, 
                           ncol = max(args$deg * estimateLowerDegrees, 1)) # out of sample average loss
  rownames(mean_test_loss) <- paste0("rep", 1:reps)
  colnames(mean_test_loss) <- paste0("deg", 1:max(args$deg * estimateLowerDegrees, 1))
  
  degToFit <- if(estimateLowerDegrees) 1:args$deg else args$deg
  
  if(args$noisy) cat("##############################\nBeginning cross-validations for degree(s)", 
                     degToFit, "and maxInteractDeg", 
                     if('maxInteractDeg' %in% names(args)) args$maxInteractDeg else args$deg,
                     "\n\n\n")
  
  
  for(i in 1:reps){
    
    test <- sample(n, nHoldout)
    ytest <- xy[test, yCol]
    args$xy <- xy[-test, ]
    
    for(j in degToFit){
      
      args$deg <- degToFit[j]
      
      trained <- do.call(polyFit, args)
      
      yhat <- predict(trained, newdata = xy[test, ])
      lost <- loss(ytest, yhat)
      if(sum(is.na(lost))){
        message("Prediction produced ", sum(is.na(lost)), 
                " NAs. Unable to obtain accurate average loss estimate for this iteration.\n\n")
      }
        
      mean_test_loss[i, j] <- mean(lost, na.rm=TRUE)
          
      if(args$noisy)
        cat("\n\nFinished iteration", i, "degree", degToFit[j], "\n", 
            "\n\t\tMean loss on test data:", mean_test_loss[i], "\n")
    }
    
  }
  
  if(length(stat) == 1){
    typical_loss <- apply(mean_test_loss, 2, stat)
  }else{
    typical_loss <- list()
    for(i in 1:length(stat))
      typical_loss[[i]] <- apply(mean_test_loss, 2, stat[[i]])
  }
  
  if(args$noisy)
    cat("\n\n\n\nFinished.\n\nTypical loss (as defined by user inputed stat) over repetitions:", 
        typical_loss, "\n\n")
  return(typical_loss)
  
}


loss <- function(y, yhat){
  
  if(is.factor(y)){
    
    if(!is.factor(yhat)){
      # assumes that characters are something like "1" for level 1,
      # not "setosa" 
      # which is the case for the predictions based on polyFit()
      # if the other type comes up, call factor(yhat, levels=levels(y))
      # before passing the vector here
      yhat <- as.numeric(yhat)
      yhat <- levels(y)[yhat]
      
    }
    
    y == yhat
  }else{
    abs(y - yhat)
  }
}  
