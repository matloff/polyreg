#################################
# FSR: Forward Stepwise Regression ###
#################################
# FSR() by default, uses 90% of data for training, 20% of which is set aside for validation.
# FSR() is primarily intended for dense matrices.
# If you have a data frame that contains factors with many levels,
# you may wish to use another function found in the polyreg library
# which takes advantage of sparse data.
# continues making the model more complicated (adding interactions, polynomials, etc.)
# until either max_poly_degree and max_interaction_degree are reached or
# improvements do not add at least threshold (default 0.01) to explained variance (of validation data).
#' FSR
#' @param Xy matrix or data.frame; outcome must be in final column.
#' @param max_poly_degree highest power to raise continuous features; default 3 (cubic).
#' @param model_type Either 'lm' (linear model), 'glm' (generalized lm), or NULL (auto-detect based on response).
#' @param max_interaction_degree highest interaction order; default 2 (allow x_i*x_j). Also interacts each level of factors with continuous features.
#' @param cor_type correlation to be used for adjusted R^2; pseudo R^2 for classification. Default "pearson"; other options "spearman" and "kendall".
#' @param threshold_include minimum improvement to include a recently added term in the model (pseudo R^2 so scale 0 to 1). -1.001 means 'include all'. Default: 0.01.
#' @param threshold_estimate minimum improvement to keep estimating (pseudo R^2 so scale 0 to 1). -1.001 means 'estimate all'. Default: 0.001.
#' @param standardize if TRUE (not default), standardizes continuous variables.
#' @param pTraining portion of data for training
#' @param pValidation portion of data for validation
#' @param min_models minimum number of models to estimate. Defaults to the number of features (unless P > N).
#' @param max_fails maximum number of models to FSR() can fail on computationally before exiting. Default == 2.
#' @param file_name If a file name (and path) is provided, saves output after each model is estimated as an .RData file. ex: file_name = "results.RData"
#' @param max_block Most of the linear algebra is done recursively in blocks to ease memory managment. Default 250. Changing up or down may slow things...
#' @param noisy display measures of fit, progress, etc. Recommended.
#' @param seed Automatically set but can also be passed as paramater.
#' @return list with slope coefficients, model details, and measures of fit
#' @export
FSR <- function(Xy,
                max_poly_degree = 3, max_interaction_degree = 2,
                model_type = NULL,
                cor_type = "pearson",
                threshold_include = 0.01,
                threshold_estimate = 0.001,
                standardize = FALSE,
                pTraining = 0.8, pValidation = 0.2,
                min_models = NULL,
                max_fails = 2,
                file_name = NULL,
                max_block = 250,
                noisy = TRUE, seed = NULL){

  if(!is.matrix(Xy) && !is.data.frame(Xy))
    stop("Xy must be a matrix or data.frame. Either way, y must be the final column.")
  if(min(pTraining, pValidation) < 0 || max(pTraining, pValidation) > 1)
    stop("pTraining and pValidation should all be between 0 and 1 and sum to 1.")
  stopifnot(is.numeric(threshold_estimate))
  stopifnot(is.numeric(threshold_include))

  out <- list()
  class(out) <- "FSR" # nested list, has S3 method
  out[["standardize"]] <- standardize
  out[["N"]] <- n <- nrow(Xy)

  Xy <- as.data.frame(Xy)
  P_factor <- 0
  factor_features <- c() # stores individual levels, omitting one

  for(i in 1:ncol(Xy)){

    k <- N_distinct(Xy[,i])

    if(k == 2 || is.character(Xy[,i])){
        Xy[,i] <- as.factor(Xy[,i])    # switch this to as.character() in case of nuissance
    }

    if(is.factor(Xy[,i])){
      P_factor <- P_factor + k - 1
      tmp <- paste(colnames(Xy)[i], "==", paste0("\'", levels(Xy[,i])[-1], "\'"))
      tmp <- paste0("(", tmp, ")")
      factor_features <- c(factor_features, tmp)
    }else{
      if(standardize)
        Xy[,i] <- scale(Xy[,i])
    }
  }

  if(is.null(model_type)){
    out[["model_type"]] <- if(is.factor(Xy[,ncol(Xy)])) "glm" else "lm"
  }else{
    if(!(model %in% c("lm", "glm")))
      stop("model must be either 'lm' or 'glm' or 'NULL' (auto-detect based on y).")
    out[["model_type"]] <- model_type
  }
  if(out$model_type == "glm"){
    if(N_distinct(Xy[,ncol(Xy)]) > 2){
      stop("multinomial outcomes not yet implemented.")
    }
  }
  out[["y_scale"]] <- if(standardize && model == "lm") sd(Xy[,ncol(Xy)]) else 1

  continuous_features <- cf <- colnames(Xy)[-ncol(Xy)][unlist(lapply(Xy[-ncol(Xy)], is_continuous))]
  P_continuous <- length(continuous_features)
  out[["continuous_features"]] <- continuous_features
  out[["factors"]] <- colnames(Xy)[-ncol(Xy)][unlist(lapply(Xy[-ncol(Xy)], is.factor))]
  out[["y"]] <- colnames(Xy)[ncol(Xy)]
  P <- P_continuous + P_factor # count without intercept
  # P does not reflect interactions or poly
  out[["P_continuous"]] <- P_continuous
  out[["P_factor"]] <- P_factor

  for(i in 2:max_poly_degree)
    continuous_features <- c(continuous_features, paste("pow(", cf, ",", i, ")"))

  features <- c(continuous_features, factor_features)
  if(max_interaction_degree > 1){
    for(i in 2:max_interaction_degree){
      features <- c(features,
                    apply(combn(c(continuous_features, factor_features), i),
                        2, isolate_interaction, i))
    }
  }


  out[["seed"]] <- if(is.null(seed)) as.integer(runif(1, 0, 10000000)) else seed
  set.seed(out$seed)
  if(noisy) message("set seed to ", out$seed, ".\n")

  out[["split"]] <- sample(c("train", "test"), n, replace=TRUE,
                  prob = c(pTraining, pValidation))

#  if(noisy) cat("N training:", length(y_train), "\nN validation:", length(y_test), "\n\n")

  # metadata to construct all possible models based on Xy and user input stored in 'add'

  models <- data.frame(features = features, stringsAsFactors = FALSE)
  models[["estimated"]]<- FALSE # will be updated based on whether successful ...
  models[["test_adjR2"]] <- NA  # computes pseudo, most important for glm. for ols, multiple definitions of R^2 equivalent including squared correlation.
                                # applies the adjustment found in adj R^2 even though out of sample ...

  if(model == "lm"){
    models[["MAPE"]] <- NA
  }else{
    models[["test_accuracy"]] <- NA
    models[["AIC"]] <- NA
    models[["BIC"]] <- NA
  }

  models[["formula"]] <- NA
  models[["accepted"]] <- FALSE
  models[["P"]] <- NA
  out[["models"]] <- models

  out[["y_name"]] <- as.character(colnames(Xy)[ncol(Xy)])
  out[["N_train"]] <- N_train <- sum(out$split == "train")
  out[["N_test"]] <- N_test <- sum(out$split == "test")
  out[["unable_to_estimate"]] <- 0

  out[["best_formula"]] <- ""
  if(is.null(min_models))
    min_models <- min(P, N_train - 1) # 'attempt all x variables additively'

  out[["max_poly_degree"]] <- max_poly_degree
  out[["max_interaction_degree"]] <- max_interaction_degree
  out[["cor_type"]] <- cor_type
  out[["threshold_include"]] <- threshold_include
  out[["threshold_estimate"]] <- 0.001
  out[["max_fails"]] <- max_fails
  out[["min_models"]] <- min_models
  out[["file_name"]] <- file_name
  out[["max_block"]] <- 250
  out[["noisy"]] <- noisy

  if(noisy) summary(out, results_overview=FALSE)

  return(FSR_estimate(out, Xy))

}

#' predict.FSR
#' @param object FSR output. Predictions will be made based on object$best_formula unless model_to_use is provided (as an integer).
#' @param newdata New Xdata.
#' @param model_to_use Integer optionally indicating a model to use if object$best_formula is not selected. Example: model_to_use = 3 will use object$models$formula[3].
#' @return y_hat (predictions using chosen model estimates).
#' @method predict FSR
#' @export
predict.FSR <- function(object, newdata, model_to_use=NULL, noisy=TRUE){

  f <- if(is.null(model_to_use)) object$best_formula else object$models$formula[model_to_use]
  beta_hat <- if(is.null(model_to_use)) object$best_coeffs else object[[mod(model_to_use)]][["est"]]
  f <- strsplit(f, "~")[[1]][2]
  f <- formula(paste("~", f))
  X_test <- model_matrix(f = f, d = newdata, noisy = noisy)
  if(object$standardize)
    for(i in 1:ncol(X_test))
      if(N_distinct(X_test[,i]) > 2)
        X_test[,i] <- scale(X_test[,i])
  return(X_test %*% beta_hat)


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

    cat("The data contains ", object$N, " observations (N_train == ",
        object$N_train, " and N_test ==", object$N_test, "), which were split using seed",
        object$seed, ". The data contains ",
        object$P_continuous,
        " continuous features, and ", object$P_factor,
        " dummy variables. Up to ", min(length(object$features), object$N_train - 1),
        " models will be estimated. (And at least ", object$min_models,
        " models will be estimated.) Each model will add a feature, which will be included in subsequent models if it explains at least an additional ",
        object$threshold_include, " of variance on [0, 1]. ", sep="")

  }

  if(results_overview){
    if(sum(object$models$estimated) == 0){

      cat("\nNo models could be estimated, likely due to (near) singularity; returning NULL. Check for highly correlated features or factors with rarely observed levels. (Increasing pTraining may help.)\n")

    }else{
      cat("\nEstimated ", sum(object$models$estimated), " models. \n")

      if(object$model_type == "lm"){
        cat("\nThe best model has Adjusted R^2:", object$best_adjR2)
        cat("\nThe best model has Mean Absolute Predicted Error:", object$best_MAPE)
      }else{
        cat("\nThe best model has pseudo R^2 (adjusted for P and N):", object$best_adjR2)
        cat("\nThe best model has out-of-sample accuracy:", object$best_test_accuracy)
      }
      cat("\n\nThe output has a data.frame out$models that contains measures of fit and information abz each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:\n\npredict(z, newdata = Xnew, model_to_use = 3) \n\n")
    }
  }

  if(!is.null(model_number)){

    m <- model_number
    cat("The added feature", ifelse(object$models$accepted[m], "WAS", "WAS NOT"), "accepted into the model.\n\n")
    cat("Adjusted R2", object$models$adjR2[m], "\n")

    if(object$model_type == "lm"){

      cat("Mean Absolute Predicted Error (MAPE)", object$models$MAPE[m], "\n")

      if(sum(object$models$accepted[1:m]) > 1){

        cat("adjusted R^2 improvement over best model so far:",
            object$models$adjR2[m] - max(object$models$adjR2[1:(m-1)], na.rm=TRUE), "\n")
        cat("MAPE improvement over best model so far:", object$models$MAPE[m] - max(object$models$MAPE[1:(m-1)], na.rm=TRUE), "\n\n\n")

      }
    }else{

      cat("(training) AIC:",  object$models$AIC[m], "\n")
      cat("(training) BIC:",  object$models$BIC[m], "\n")
      cat("(test) classification accuracy:", object$models$test_accuracy[m], "\n")

      if(sum(object$models$accepted) > 1){
        cat("pseudo-R^2 (adjusted based on P and N) improvement over best model so far:",
            object$models$adjR2[m] - max(object$models$adjR2[1:(m-1)], na.rm=TRUE), "\n")
        cat("Classification accuracy improvement on the test data:",
            object$models$test_accuracy[m] - max(object$models$accuracy[1:(m-1)], na.rm=TRUE), "\n")
      }
    }

  }
<<<<<<< HEAD
=======
}

#' predict.FSR
#' @param object FSR output. Predictions will be made based on object$best_formula unless model_to_use is provided (as an integer).
#' @param newdata New Xdata.
#' @param model_to_use Integer optionally indicating a model to use if object$best_formula is not selected. Example: model_to_use = 3 will use object$models$formula[3].
#' @return y_hat (predictions using chosen model estimates).
#' @method predict FSR
#' @export
predict.FSR <- function(object, newdata, model_to_use=NULL, noisy=TRUE){

  f <- if(is.null(model_to_use)) object$best_formula else object$models$formula[model_to_use]
  beta_hat <- if(is.null(model_to_use)) object$best_coeffs else object[[mod(model_to_use)]][["est"]]
  f <- strsplit(f, "~")[[1]][2]
  f <- formula(paste("~", f))
  X_test <- model_matrix(f = f, d = newdata, noisy = noisy)
  if(object$standardize)
    for(i in 1:ncol(X_test))
      if(N_distinct(X_test[,i]) > 2)
        X_test[,i] <- scale(X_test[,i])
  return(X_test %*% beta_hat)
>>>>>>> 47c8b1c1581439783c67fb455337a3c184467b21


}
