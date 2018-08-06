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
#' @param max_interaction_degree highest interaction order; default 1. Also interacts each level of factors with continuous features.
#' @param cor_type correlation to be used for adjusted R^2; pseudo R^2 for classification. Default "pearson"; other options "spearman" and "kendall".
#' @param threshold_include minimum improvement to include a recently added term in the model (pseudo R^2 so scale 0 to 1). -1.001 means 'include all'. Default: 0.01.
#' @param threshold_estimate minimum improvement to keep estimating (pseudo R^2 so scale 0 to 1). -1.001 means 'estimate all'. Default: 0.001.
#' @param standardize if TRUE (not default), standardizes continuous variables.
#' @param pTraining portion of data for training
#' @param pValidation portion of data for validation
#' @param min_models minimum number of models to estimate. Defaults to the number of features.
#' @param file_name If a file name (and path) is provided, saves output after each model is estimated as an .RData file. ex: file_name = "results.RData"
#' @param max_block_size Most of the linear algebra is done recursively in blocks to ease memory managment. Default 250. Changing up or down may slow things...
#' @param noisy display measures of fit, progress, etc. Recommended.
#' @param seed Automatically set but can also be passed as paramater.
#' @return list with slope coefficients, model details, and measures of fit
#' @export
FSR <- function(Xy,
                max_poly_degree = 3, max_interaction_degree = 1,
                cor_type = "pearson",
                threshold_include = 0.01,
                threshold_estimate = 0.001,
                standardize = FALSE,
                pTraining = 0.8, pValidation = 0.2,
                min_models = NULL,
                file_name = NULL,
                max_block_size = 250,
                noisy = TRUE, seed = NULL,
                model = "lm"){

  if(!is.matrix(Xy) && !is.data.frame(Xy))
    stop("Xy must be a matrix or data.frame. Either way, y must be the final column.")
  if(min(pTraining, pValidation) < 0 || max(pTraining, pValidation) > 1)
    stop("pTraining and pValidation should all be between 0 and 1 and sum to 1.")
  stopifnot(is.numeric(threshold_estimate))
  stopifnot(is.numeric(threshold_include))

  out <- list()
  class(out) <- "FSR" # nested list, has S3 method
  out[["n"]] <- n <- nrow(Xy)

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

  if(is.null(model)){
    model <- if(is.factor(Xy[,ncol(Xy)])) "glm" else "lm"
  }else{
    if(!(model %in% c("lm", "glm")))
      stop("model must be either 'lm' or 'glm' or 'NULL' (auto-detect based on y).")
  }
  if(model == "glm"){
    if(N_distinct(Xy[,ncol(Xy)]) > 2){
      stop("multinomial outcomes not implemented.")
    }else{
      warning("Logistic regression not yet implemented, treating as continuous...")
    }
  }

  continuous_features <- cf <- colnames(Xy)[-ncol(Xy)][unlist(lapply(Xy[-ncol(Xy)], is_continuous))]
  P_continuous <- length(continuous_features)
  out[["continuous_features"]] <- continuous_features
  out[["factors"]] <- colnames(Xy)[-ncol(Xy)][unlist(lapply(Xy[-ncol(Xy)], is.factor))]
  out[["y"]] <- colnames(Xy)[ncol(Xy)]
  P <- P_continuous + P_factor # count without intercept
  # P does not reflect interactions or poly

  for(i in 2:max_poly_degree)
    continuous_features <- c(continuous_features, paste("pow(", cf, ",", i, ")"))

  features <- c(continuous_features, factor_features)
  if(max_interaction_degree > 0){
    for(i in 1:max_interaction_degree){
      features <- c(features,
                    apply(combn(c(continuous_features, factor_features), i + 1),
                        2, isolate_interaction, i))
    }
  }


  out[["seed"]] <- if(is.null(seed)) as.integer(runif(1, 0, 10000000)) else seed
  set.seed(out$seed)
  if(noisy) message("set seed to ", out$seed, ".\n")

  out[["split"]] <- sample(c("train", "test"), n, replace=TRUE,
                  prob = c(pTraining, pValidation))
  y_train <- Xy[out$split == "train", ncol(Xy)]
  y_test <- Xy[out$split == "test", ncol(Xy)]

  if(noisy) cat("N training:", length(y_train), "\nN validation:", length(y_test), "\n\n")

  # metadata to construct all possible models based on Xy and user input stored in 'add'

  models <- data.frame(added_feature = features, stringsAsFactors = FALSE)
  models[["estimated"]]<- FALSE # will be updated based on whether successful ...
  models[["adjR2"]] <- NA
  models[["MAPE"]] <- NA
  models[["formula"]] <- NA
  models[["accepted"]] <- FALSE
  out[["models"]] <- models

  m <- 1            # counts which model
  improvement <- 0  # not meaningful; just initializing ...
  y_name <- as.character(colnames(Xy)[ncol(Xy)])
  out[["N_train"]] <- N_train <- sum(out$split == "train")
  out[["N_test"]] <- N_test <- sum(out$split == "test")
  unable_to_estimate <- 0 # allowed 2 fails ... change?
  out[["best_formula"]] <- ""
  if(is.null(min_models))
    min_models <- P # 'attempt all x variables additively'

  if(noisy) cat("The data contains", n, "observations (N_train ==", N_train, "and N_test ==", N_test, "),",
                P_continuous,
                "continuous features, and", P_factor,
                "dummy variables. Up to", length(features),
                "models will be estimated. (And at least", min_models, "models will be estimated.) Each model will add a feature, which will be included in subsequent models if it explains at least an additional", threshold_include, "on the adjusted R^2 scale.")



  while((m <= nrow(models)) &&
        ((improvement > threshold_estimate) || m <= min_models) &&
        unable_to_estimate < 2){

    out[[mod(m)]] <- list()

    if(sum(out$models$accepted) == 0){
      out$models$formula[m] <- paste(y_name, "~", features[m])
    }else{
      out$models$formula[m] <- paste(out[["best_formula"]], "+", features[m])
    }
    X_train <- model_matrix(formula(out$models$formula[m]), Xy[out$split == "train",], noisy=TRUE)

    if(!exists("X_train")){

      if(noisy) message("Unable to construct model.matrix for model ",  m, ". Skipping.")
      unable_to_estimate <- unable_to_estimate + 1

    }else{

      out[[mod(m)]][["p"]] <- ncol(X_train)

      if(out[[mod(m)]][["p"]] >= nrow(X_train)){ # N_train

        if(noisy) message("There are too few training observations to estimate model ",  m, ". Skipping.")
        unable_to_estimate <- unable_to_estimate + 1

      }else{

        if(sum(out$models$accepted) == 0){
          XtX_inv <- block_solve(X = X_train,  max_block_size = max_block_size)
        }else{
          XtX_inv <- block_solve(X = X_train, A_inv = XtX_inv_accepted, max_block_size = max_block_size)
        }
        # passing X takes crossproduct first
        # starting with second iteration,
        # XtX_inv of the last accepted model is taken as the inverse of the first block


        if(!is.null(XtX_inv)){

          out[[mod(m)]][["est"]] <- tcrossprod(XtX_inv, X_train) %*% y_train

          if(m == length(features))
            remove(XtX_inv)

          if(sum(is.na(out[[mod(m)]][["est"]])) > 0){

            if(noisy) cat("\nfailed to estimate model", m, "skipping...\n")
            unable_to_estimate <- unable_to_estimate + 1

          }else{

            out$models$estimated[m] <- TRUE
            out[[mod(m)]][["y_hat"]] <- model_matrix(formula(out$models$formula[m]),
                                                     Xy[out$split == "test", ]) %*% out[[mod(m)]][["est"]]
            R2 <- cor(out[[mod(m)]][["y_hat"]], y_test, method=cor_type)^2
            adjR2 <- (N_train - ncol(X_train) - 1)/(N_train - 1)*R2
            out$models$adjR2[m] <- out[[mod(m)]][[paste0("adj_R2_", cor_type)]] <- adjR2

            MAPE <- mean(abs(out[[mod(m)]][["y_hat"]] - y_test))
            out$models$MAPE[m] <- out[[mod(m)]][["MAPE"]] <- MAPE

            if(sum(out$models$accepted) == 0){
              improvement <- adjR2
            }else{
              improvement <- adjR2 - out$best_adjR2
            }

            if(improvement > threshold_include){


              out$best_formula <- out$models$formula[m]
              out[["best_coeffs"]] <- out[[mod(m)]][["est"]]
              out[["best_adjR2"]] <- adjR2
              out[["best_MAPE"]] <- MAPE
              out$models$accepted[m] <- TRUE
              XtX_inv_accepted <- XtX_inv

            }

            if(noisy){
              cat("\n\n\n\n")
              cat("Model:", m, "\n\n")
              cat("The added feature", ifelse(out$models$accepted[m], "WAS", "WAS NOT"), "accepted into the model.")
              cat("Adjusted R2", out$models$adjR2[m], "\n")
              cat("Mean Absolute Predicted Error (MAPE)", out$models$MAPE[m], "\n")
              if(sum(out$models$accepted) > 0){
                cat("adj R2 improvement over best model so far:", improvement, "\n")
                cat("MAPE improvement over best model so far:", MAPE - out$best_MAPE, "\n\n\n")
              }
            }
          }
        }else{

          if(noisy) cat("\nunable to estimate Model", m, "likely due to (near) singularity.\n")
          unable_to_estimate <- unable_to_estimate + 1

        }
      }

    }

    m <- m + 1
    if(!is.null(file_name))
      save(out, file=file_name)

  } # end WHILE loop


  if(sum(out$models$estimated) == 0){

    if(noisy) message("\nNo models could be estimated, likely due to (near) singularity; returning NULL. Check for highly correlated features or factors with rarely observed levels. (Increasing pTraining may help.)\n")
    return(NULL)

  }else{
    if(noisy) message("\nEstimated ", sum(out$models$estimated), " models. The output will have a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list object, such that out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:\n\npredict(out, newdata = Xnew, model_to_use = 3) \n\n\n")
    return(out)
  }
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
  return(X_test %*% beta_hat)


}


