
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

rep_kms <- function(xy, input_formula, reps = 5, stat = median, 
                    nHoldout = min(10000, round(0.2*nrow(xy))),
                    yCol = NULL,
                    seed = list(seed = NULL, disable_gpu = FALSE, disable_parallel_cpu = FALSE), 
                    noisy = FALSE, ...){
  
  require(kerasformula)
  
  args <- list(...)
  args$verbose <- noisy
  args$input_formula <- input_formula
  
  # handle using predict() instead, override user preferences here...
  args$pTraining <- 1
  args$validation_split <- 0
  
  if(!is.null(yCol)){
    xy <- cbind(xy[,-yCol], xy[,yCol])
  }else{
    yCol <- ncol(xy)
  }
  
  if(!is.null(seed$seed)){
    args$seed <- seed
    set.seed(seed$seed)
    if(noisy) 
      message("seed set to ", seed$seed, " using ", sessionInfo()[[1]]$version.string, 
      ".
      
NOTE: R upgraded its integer sampling method as version R 3.6.0, so results may differ slightly depending...

NOTE: keras calls languages external to R. 
To fully remove simulation error set, for example,

list(seed = 1234, disable_gpu = TRUE, disable_parallel_cpu = TRUE). ")
  }
  
  n <- nrow(xy)
  
  mean_test_loss <- matrix(nrow = reps, 
                           ncol = 1) # out of sample average loss
  rownames(mean_test_loss) <- paste0("rep", 1:reps)
  
  for(i in 1:reps){
    
    test <- sample(n, nHoldout)
    args$data <- xy[-test,]
    trained <- do.call(kms, args)
    predicted <- predict(trained, newdata = xy[test, ])
    mean_test_loss[i] <- mean(loss(predicted$y_test, predicted$fit))
    cat("\n\nmean absolute predicted error on iteration", i, "is",  mean_test_loss[i])
  }
  if(length(stat) == 1){
    return(stat(mean_test_loss))
  }else{
    out <- list()
    for(i in 1:length(stat))
      out[[i]] <- stat[[i]](mean_test_loss)
    return(out)
  }
}

rep_xvalPoly <- function(xy, deg, reps = 5, stat = median, 
                         nHoldout = min(10000, round(0.2*nrow(xy))), 
                         yCol = NULL, 
                         seed = NULL, estimateLowerDegrees = TRUE, ...){
  
  args <- list(...)
  args$deg <- deg
  args$noisy <- if(is.null(args$noisy)) TRUE else args$noisy
  
  
  n <- nrow(xy)
  
  if(!is.null(yCol)){
    xy <- cbind(xy[,-yCol], xy[,yCol])
  }else{
    yCol <- ncol(xy)
  }
    
  if(length(unique(xy[,yCol])) == 2 || is.character(xy[,yCol]))
    xy[,yCol] <- as.factor(xy[,yCol])
  
  if(!is.null(seed)){
    set.seed(seed)
    if(noisy) message("seed set to ", seed, " using ", sessionInfo()[[1]]$version.string, 
                      ".
                      
NOTE: R upgraded its integer sampling method as version R 3.6.0, so results may differ slightly depending...")
  }
  
  
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
