#' @export
predict.FSR <- function(object, newdata, model_to_use=NULL, 
                        standardize=NULL, noisy=TRUE,...)
{
  

  if(!is.null(standardize) && object$standardize){
    for(var_name in names(object[["train_scales"]])){
      if(var_name %in% colnames(newdata)){
        newdata[[var_name]] <- (newdata[[var_name]] - object[["train_scales"]][[var_name]][["mean"]])/object[["train_scales"]][[var_name]][["sd"]]
      }
    }
  }

  m <- if(is.null(model_to_use)) which(object$best_formula == object$models$formula) else model_to_use

  mf <- if(is.null(m)) object$best_formula else object$models$formula[m]

  mf <- strsplit(mf, "~")[[1]][2]
  mf <- formula(paste("~", mf))
  X_test <- model_matrix(modelFormula = mf, dataFrame = newdata, noisy = noisy, intercept = TRUE)

  y_hat <- X_test %*% object[[mod(m)]][["coeffs"]]

  if(object$outcome == "continuous"){

    return(y_hat) # should this be coeffs?

  }else{

    if(object$outcome == "binary"){

      return(list(probs = y_hat,
                  classified = classify(y_hat,
                                        labels = object$training_labels,
                                        cutoff = object$y_train_mean)))

    }else{  # multinomial case

      p_K <- 1/(1 + rowSums(exp(y_hat))) # prob of being in reference category
      pred <- list(probs = cbind(p_K, exp(y_hat)*p_K))
      colnames(pred$probs)[1] <- object$reference_category
      pred[["classified"]] <- classify(pred$probs)

    }

    return(pred)
  }

}

#' @export
print.polyFit <- function(x, ...){
  print(x$fit)
}

#' @export
print.FSR <- function(x, ...){
  summary(x, estimation_overview=FALSE, results_overview=TRUE)
}

#' @export
summary.FSR <- function(object, estimation_overview=TRUE, results_overview=TRUE, model_number = NULL, ...){

  if(estimation_overview){

    message("The dependent variable is '",  object$y_name,
        "' which will be treated as ",  object$outcome, ". ",
        " The data contains ", object$N, " observations (N_train == ",
        object$N_train, " and N_test == ", object$N_test,
        "), which were split using seed ", object$seed, ". The data contains ",
        object$P_continuous,
        " continuous features and ", object$P_factor,
        " dummy variables. Between ", object$min_models, " and ",
        min(nrow(object$models), object$N_train - 1),
        " models will be considered. Each model will add a feature, ",
        "which will be included in subsequent models if it explains at least an additional ",
        object$threshold_include,
        " of variance out-of-sample (after adjusting for the additional term on [0, 1]).",
        " Note: x3*x7, for example, will not be included if x3 and and x7 were not previously included.\n\n",
        sep="")

  }

  if(results_overview){

    if(sum(object$models$estimated) == 0){

      message("\nNo models could be estimated, likely due to (near) singularity; returning NULL. Check for highly correlated features or factors with rarely observed levels. (Increasing pTraining may help.)\n")

    }else{
      message("\nEstimated ", sum(object$models$estimated), " models. \n")

      if(object$outcome == "continuous"){
        message("\nThe best model has Predicted Adjusted R^2: ", round(object$best_adjR2, 4))
        message("\nThe best model has Mean Absolute Predicted Error: ", round(object$best_MAPE, 4))
      }else{
        if(object$outcome == "binary"){
          message("\nThe best model has pseudo R^2 (adjusted for P and N): ",
              object$best_test_adjR2)
          message("\nThe best model has out-of-sample accuracy: ",
              object$models$test_accuracy[which(object$best_formula == object$models$formula)])
        }else{
          message("\nThe best model has out-of-sample accuracy (adjusted for P and N):",
              object$best_test_adj_accuracy)
        }
      }
      message("\n\nThe output has a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:\n\npredict(z, newdata = Xnew, model_to_use = 3) \n\n")
    }
  }

  if(!is.null(model_number)){

    m <- model_number
    message("\n\n\nThe added coefficient, corresponding to ", object$models$features[m],
        ifelse(object$models$accepted[m], ", WAS", ", WAS NOT"),
        " accepted into model ", m, ".\n\n", sep="")

    if(object$outcome == "continuous"){

      message("Adjusted R^2: ", round(object$models$test_adjR2[m], 5), "\n")
      message("Mean Absolute Predicted Error (MAPE): ", round(object$models$MAPE[m], 5), "\n")

      if(sum(object$models$accepted[1:m]) > 1){

        message("Adjusted R^2 improvement over best model so far: ",
            round(object$models$test_adjR2[m] - max(object$models$test_adjR2[1:(m - 1)], na.rm=TRUE), 5),
            "\n")
        message("MAPE improvement over best model so far: ",
            round(object$models$MAPE[m] - max(object$models$MAPE[1:(m - 1)], na.rm=TRUE), 5),
            "\n\n\n")

      }
    }else{

      if(object$outcome == "binary")
        message("Pseudo R^2 (adjusted for P and N): ", round(object$models$test_adjR2[m], 5), "\n")

      if(!object$linear_estimation){
        message("(training) AIC: ",  round(object$models$AIC[m], 5), "\n")
        message("(training) BIC: ",  round(object$models$BIC[m], 5), "\n")
      }
      message("(test) classification accuracy:", round(object$models$test_accuracy[m], 5), "\n")

      if(sum(object$models$accepted) > 1){
        message("Classification accuracy improvement on the test data: ",
            round(object$models$test_accuracy[m] - max(object$models$test_accuracy[1:(m - 1)], na.rm=TRUE), 5), 
            "\n")
        if(object$outcome == "binary"){
          message("pseudo-R^2 (adjusted based on P and N) improvement over best model so far:",
              round(object$models$test_adjR2[m] - max(object$models$test_adjR2[1:(m - 1)], na.rm=TRUE), 5),
              "\n")
        }
      }
    }
  }
}

##################################################################
# predict.polyFit: predict the fitted models on newdata
##################################################################

predict.polyFit <- function(object, newdata, ...)
{
  use <- object$use
  
  # the next couple dozen lines are devoted to forming plm.newdata, which
  # will ultimately be fed into predict.lm(), predict.glm() or whatever;
  # to do this, newdata, the argument above, must be expanded to
  # polynomial form, and possibly run through PCA
  
  doPCA <- !is.null(object$PCA)
  
  if (!doPCA) {
    
    plm.newdata <- getPoly(newdata, object$degree, 
                           object$maxInteractDeg, 
                           modelFormula = object$XtestFormula,
                           retainedNames = object$retainedNames,
                           ...)$xdata
    
  } else if (object$PCA == "prcomp") {
    
    message("Beginning PCA\n\n", timestamp())
    
    if (object$pcaLocation == "front") {
      
      new_data <- predict(object$pca.xy, newdata)[,1:object$pcaCol]
      plm.newdata <- getPoly(new_data, object$degree, object$maxInteractDeg, ...)$xdata
      
    } else if (object$pcaLocation == "back") {
      
      new_data <- getPoly(newdata, object$degree, object$maxInteractDeg, ...)$xdata
      plm.newdata <- predict(object$pca.xy, new_data)[,1:object$pcaCol]
      plm.newdata <- as.data.frame(plm.newdata)
      
    } else stop('invalid pcaLocation. Should be "front" or "back".')
  } else if (object$PCA == "RSpectra") {
    
    if (object$pcaLocation == "front") {
      
      xy.eig <- object$pca.xy
      new_data <- as.matrix(newdata) %*% xy.eig$vectors
      plm.newdata <- getPoly(new_data, object$degree, object$maxInteractDeg)$xdata
      
    } else if (object$pcaLocation == "back") {
      new_data <- getPoly(newdata, object$degree, object$maxInteractDeg, ...)$xdata
      ### xy.cov <- cov(new_data)
      ### xy.eig <- eigs(xy.cov,object$pcaCol)
      xy.eig <- object$pca.xy
      plm.newdata <-
        as.matrix(new_data) %*% xy.eig$vectors[,1:object$pcaCol]
      plm.newdata <- as.data.frame(plm.newdata)
    }
  message("Finished with PCA and model matrix construction.\n\n", timestamp())
  }  # end doPCA

  
  
  if (object$use == "lm") {
    pred <- predict(object$fit, plm.newdata)
    return(pred)
  }
  
  if (object$use == "mvrlm") {
    pre <- predict(object$fit, plm.newdata)
    pred <- apply(pre,1,which.max)
    return(pred)
  }
  
  # glm case
  if (is.null(object$glmMethod)) { # only two classes
    pre <- predict(object$fit, plm.newdata, type = 'response')
    pred <- ifelse(pre > 0.5, object$classes[1], object$classes[2])
  } else { # more than two classes
    len <- length(object$classes)
    if (object$glmMethod == "multlog") { # multinomial logistics
      pr <- predict(object$fit, plm.newdata, type="probs")
      idx <- apply(pr,1, which.max)
      col.name <- colnames(pr)
      lc <- length(col.name)
      tempM <- matrix(rep(col.name, length(idx)), ncol=lc, byrow = TRUE)
      pred <- NULL
      for (r in 1:nrow(tempM)) {
        pred[r] <- tempM[r,idx[r]]
      }
      return(pred)
    } # end multinomial logistics
    else if (object$glmMethod == "all") { # all-vs-all method
      votes <- matrix(0, nrow = nrow(plm.newdata), ncol = len)
      for (i in 1:len) {
        for (j in 1:len) {
          if (i == j)
            next
          pre <- predict(object$fit[[i]][[j]], plm.newdata, type="response")
          votes[,i] <- votes[,i] + ifelse(pre > 0.5, 1, 0)
        } # for i
      } # for j
      winner <- apply(votes, 1, which.max)
      
    } else if (object$glmMethod == "one") { # one-vs-all method
      prob <- matrix(0, nrow=nrow(plm.newdata), ncol=len)
      for (i in 1:len) {
        # prob[,i] <- parSapplyLB(object$fit[[i]],
        prob[,i] <- predict(object$fit[[i]],
                            plm.newdata, type = "response")
      }
      winner <- apply(prob, 1, which.max)
    } # one-vs-all method
    # calculate pred for all-vs-all & one-vs-all
    pred <- NULL
    for (k in 1:nrow(plm.newdata)) {
      pred[k] <- object$classes[winner[k]]
    }
  } # end more than two classes
  return(pred)
  # end glm case
  
}
