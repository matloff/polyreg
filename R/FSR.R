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
#' @param threshold minimum improvement to keep estimating (pseudo R^2 so scale 0 to 1). -1.001 means 'estimate all'. Default: 0.01.
#' @param standardize if TRUE (default), standardizes continuous variables.
#' @param pTraining portion of data for training
#' @param pValidation portion of data for validation
#' @param max_block_size Most of the linear algebra is done recursively in blocks to ease memory managment. Default 250. Changing up or down may slow things...
#' @param noisy display measures of fit, progress, etc. Recommended.
#' @param seed Automatically set but can also be passed as paramater.
#' @return list with slope coefficients, model details, and measures of fit
#' @export
FSR <- function(Xy,
                max_poly_degree = 3, max_interaction_degree = 1,
                cor_type = "pearson", threshold = 0.01, standardize = TRUE,
                pTraining = 0.8, pValidation = 0.2, max_block_size = 250,
                noisy = TRUE, seed = NULL,
                model = "lm"){

  if(!is.matrix(Xy) && !is.data.frame(Xy))
    stop("Xy must be a matrix or data.frame. Either way, y must be the final column.")
  if(min(pTraining, pValidation) < 0 || max(pTraining, pValidation) > 1)
    stop("pTraining and pValidation should all be between 0 and 1 and sum to 1.")
  stopifnot(is.numeric(threshold))

  out <- list()
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

  if(noisy) cat("The data contains", n, "observations,", length(continuous_features),
                "continuous features, and", P_factor,
                "dummy variables. Up to", length(features),
                "models will be estimated. The output will have a data.frame that contains measures of fit and information about each model, such as the formula call. The output is also a nested list object, such that out$model1, out$model2, etc. contains further metadata.\n\n")


  m <- 1            # counts which model
  improvement <- 0  # not meaningful; just initializing ...
  y_name <- as.character(colnames(Xy)[ncol(Xy)])
  out[["N_train"]] <- N_train <- sum(out$split == "train")
  out[["N_test"]] <- N_test <- sum(out$split == "train")
  unable_to_estimate <- 0 # allowed 2 fails ... change?
  out[["best_formula"]] <- ""

  while((m <= nrow(models)) &&
        ((improvement > threshold) || m <= P) && # 'attempt all x variables additively'
        unable_to_estimate < 2){

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

      out[[mod(m)]] <- list()

      out[[mod(m)]][["p"]] <- ncol(X_train)

      if(out[[mod(m)]][["p"]] >= nrow(X_train)){

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

          out[[mod(m)]][["y_hat"]] <- model_matrix(formula(out$models$formula[m]),
                                                   Xy[out$split == "test", ]) %*% out[[mod(m)]][["est"]]
          R2 <- cor(out[[mod(m)]][["y_hat"]], y_test, method=cor_type)^2
          out[[mod(m)]][[paste0("adj_R2_", cor_type)]] <- (length(y_train) - ncol(X_train) - 1)/(length(y_train) - 1)*R2
          out[[mod(m)]][["MAPE"]] <- mean(abs(out[[mod(m)]][["y_hat"]] - y_test))

          improvement <- if(m == 1) out[[mod(m)]][[paste0("adj_R2_", cor_type)]] else out[[mod(m)]][[paste0("adj_R2_", cor_type)]] - out[[mod(m - 1)]][[paste0("adj_R2_", cor_type)]]

          out$models$estimated[m] <- TRUE
          out$models$adjR2[m] <- out[[mod(m)]][[paste0("adj_R2_", cor_type)]]
          out$models$MAPE[m] <- out[[mod(m)]][["MAPE"]]

          if(noisy){
            cat("\n\n")
            cat("model", m, "adjusted R2", out$models$adjR2[m], "\n")
            cat("model", m, "Mean Absolute Predicted Error (MAPE)", out$models$MAPE[m], "\n")
            if(m > 1){
              cat("adj R2 improvement over last model:", out$models$adjR2[m] - out$models$adjR2[m - 1], "\n")
              cat("MAPE improvement over last model:", out$models$MAPE[m] - out$models$MAPE[m - 1], "\n")
            }

            cat("\n\n")
          }

          if(m > 1 && (out$models$adjR2[m] > out$models$adjR2[m - 1])){
            out$best_formula <- out$models$formula[m]
            out$models$accepted[m] <- TRUE
          }else{

            if(out$models$adjR2[m] > threshold){
              out$best_formula <- out$models$formula[m]
              out$models$accepted[m] <- TRUE
            }

          }
          if(out$models$accepted[m])
            XtX_inv_accepted <- XtX_inv


        }else{

          if(noisy) cat("\nunable to estimate Model", m, "likely due to (near) singularity.\n")
          unable_to_estimate <- unable_to_estimate + 1

        }
      }

    }

    m <- m + 1

  } # end WHILE loop


  if(sum(out$models$estimated) == 0){

    if(noisy) message("\nNo models could be estimated, likely due to (near) singularity; returning NULL. Check for highly correlated features or factors with rarely observed levels. (Increasing pTraining may help.)\n")
    return(NULL)

  }else{
    return(out)
  }
}
