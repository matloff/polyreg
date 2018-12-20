#' predict.FSR
#' @param object FSR output. Predictions will be made based on object$best_formula unless model_to_use is provided (as an integer).
#' @param newdata New Xdata.
#' @param model_to_use Integer optionally indicating a model to use if object$best_formula is not selected. Example: model_to_use = 3 will use object$models$formula[3].
#' @param standardize Logical--standardize numeric variables? (If NULL, the default, bypasses and decides based on object$standardize.)
#' @param noisy Display output?
#' @return y_hat (predictions using chosen model estimates).
#' @method predict FSR
#' @examples
#' out <- FSR(mtcars[1:20,])
#' forecast <- predict(out, mtcars[21:nrow(mtcars),])
#' @export
predict.FSR <- function(object, newdata, model_to_use=NULL, standardize=NULL, noisy=TRUE){


  if(!is.null(standardize) && object$standardize){
    for(var_name in names(object[["train_scales"]])){
      if(var_name %in% colnames(newdata)){
        newdata[[var_name]] <- (newdata[[var_name]] - object[["train_scales"]][[var_name]][["mean"]])/object[["train_scales"]][[var_name]][["sd"]]
      }
    }
  }

  m <- if(is.null(model_to_use)) which(object$best_formula == object$models$formula) else model_to_use

  f <- if(is.null(m)) object$best_formula else object$models$formula[m]

  f <- strsplit(f, "~")[[1]][2]
  f <- formula(paste("~", f))
  X_test <- model_matrix(f = f, d = newdata, noisy = noisy, intercept = TRUE)

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


#' summary.FSR
#' @param object an FSR object
#' @param estimation_overview logical: describe how many models were planned, sample size, etc.?
#' @param results_overview logical: give overview of best fit model, etc?
#' @param model_number If non-null, an integer indicating which model to display a summary of.
#' @method summary FSR
#' @export
summary.FSR <- function(object, estimation_overview=TRUE, results_overview=TRUE, model_number = NULL){

  if(estimation_overview){

    cat("The dependent variable is '",  object$y_name,
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

      cat("\nNo models could be estimated, likely due to (near) singularity; returning NULL. Check for highly correlated features or factors with rarely observed levels. (Increasing pTraining may help.)\n")

    }else{
      cat("\nEstimated ", sum(object$models$estimated), " models. \n")

      if(object$outcome == "continuous"){
        cat("\nThe best model has Predicted Adjusted R^2:", object$best_adjR2)
        cat("\nThe best model has Mean Absolute Predicted Error:", object$best_MAPE)
      }else{
        if(object$outcome == "binary"){
          cat("\nThe best model has pseudo R^2 (adjusted for P and N):",
              object$best_test_adjR2)
          cat("\nThe best model has out-of-sample accuracy:",
              object$models$test_accuracy[which(object$best_formula == object$models$formula)])
        }else{
          cat("\nThe best model has out-of-sample accuracy (adjusted for P and N):",
              object$best_test_adj_accuracy)
        }
      }
      cat("\n\nThe output has a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:\n\npredict(z, newdata = Xnew, model_to_use = 3) \n\n")
    }
  }

  if(!is.null(model_number)){

    m <- model_number
    cat("\n\n\nThe added coefficient, corresponding to ", object$models$features[m],
        ifelse(object$models$accepted[m], ", WAS", ", WAS NOT"),
        " accepted into model ", m, ".\n\n", sep="")

    if(object$outcome == "continuous"){

      cat("Adjusted R2", object$models$test_adjR2[m], "\n")
      cat("Mean Absolute Predicted Error (MAPE)", object$models$MAPE[m], "\n")

      if(sum(object$models$accepted[1:m]) > 1){

        cat("adjusted R^2 improvement over best model so far:",
            object$models$test_adjR2[m] - max(object$models$test_adjR2[1:(m - 1)], na.rm=TRUE),
            "\n")
        cat("MAPE improvement over best model so far:",
            object$models$MAPE[m] - max(object$models$MAPE[1:(m - 1)], na.rm=TRUE),
            "\n\n\n")

      }
    }else{

      if(object$outcome == "binary")
        cat("Pseudo R2 (adjusted for P and N)", object$models$test_adjR2[m], "\n")

      if(!object$linear_estimation){
        cat("(training) AIC:",  object$models$AIC[m], "\n")
        cat("(training) BIC:",  object$models$BIC[m], "\n")
      }
      cat("(test) classification accuracy:", object$models$test_accuracy[m], "\n")

      if(sum(object$models$accepted) > 1){
        cat("Classification accuracy improvement on the test data:",
            object$models$test_accuracy[m] - max(object$models$test_accuracy[1:(m - 1)], na.rm=TRUE),
            "\n")
        if(object$outcome == "binary"){
          cat("pseudo-R^2 (adjusted based on P and N) improvement over best model so far:",
              object$models$test_adjR2[m] - max(object$models$test_adjR2[1:(m - 1)], na.rm=TRUE),
              "\n")
        }
      }
    }
  }
}
