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
#' @param outcome Treat y as either 'continuous', 'binary', 'multinomial', or NULL (auto-detect based on response).
#' @param linear_estimation Logical: model outcome as linear and estimate with ordinary least squares? Recommended for speed on large datasets even if outcome is categorical. (For multinomial outcome, this means treated response as vector.) If FALSE, estimator chosen based on 'outcome' (i.e., OLS for continuous outcomes, glm() to estimate logistic regression models for 'binary' outcomes, and nnet::multinom() for 'multinomial').
#' @param max_interaction_degree highest interaction order; default 2 (allow x_i*x_j). Also interacts each level of factors with continuous features.
#' @param threshold_include minimum improvement to include a recently added term in the model (change in fit originally on 0 to 1 scale). -1.001 means 'include all'. Default: 0.01. (Adjust R^2 for linear models, Pseudo R^2 for logistic regression, out-of-sample accuracy for multinomial models. In latter two cases, the same adjustment for number of predictors is applied as pseudo-R^2.)
#' @param threshold_estimate minimum improvement to keep estimating (pseudo R^2 so scale 0 to 1). -1.001 means 'estimate all'. Default: 0.001.
#' @param min_models minimum number of models to estimate. Defaults to the number of features (unless P > N).
#' @param max_fails maximum number of models to FSR() can fail on computationally before exiting. Default == 2.
#' @param standardize if TRUE (not default), standardizes continuous variables.
#' @param pTraining portion of data for training
#' @param file_name If a file name (and path) is provided, saves output after each model is estimated as an .RData file. ex: file_name = "results.RData". See also store_fit for options as to how much to store in the outputted object.
#' @param store_fit If file_name is provided, FSR() will return coefficients, measures of fit, and call details. Save entire fit objects? Options include "none" (default, just save those other items), "accepted_only" (only models that meet the threshold), and "all".
#' @param max_block Most of the linear algebra is done recursively in blocks to ease memory managment. Default 250. Changing up or down may slow things...
#' @param noisy display measures of fit, progress, etc. Recommended.
#' @param seed Automatically set but can also be passed as paramater.
#' @return list with slope coefficients, model details, and measures of fit
#' @examples
#' out <- FSR(mtcars)
#' @importFrom nnet multinom
#' @importFrom stats relevel glm
#' @export
FSR <- function(Xy,
                max_poly_degree = 3, max_interaction_degree = 2,
                outcome = NULL, linear_estimation = FALSE,
                threshold_include = 0.01, threshold_estimate = 0.001,
                min_models = NULL, max_fails = 2,
                standardize = FALSE,
                pTraining = 0.8,
                file_name = NULL,
                store_fit = "none",
                max_block = 250,
                noisy = TRUE, seed = NULL){

  if(!is.matrix(Xy) && !is.data.frame(Xy))
    stop("Xy must be a matrix or data.frame. Either way, y must be the final column.")
  if(pTraining <= 0 || pTraining > 1)
    stop("pTraining should all be between 0 and 1.")
  pValidation <- 1 - pTraining
  stopifnot(is.numeric(threshold_estimate))
  stopifnot(is.numeric(threshold_include))
  stopifnot(is.logical(linear_estimation))

  store_fit <- match_arg(store_fit, c("none", "accepted", "all"))
  outcome <- match_arg(outcome, c("continuous", "binary", "multinomial"))

  out <- list()
  class(out) <- "FSR" # nested list, has S3 method
  out[["standardize"]] <- standardize
  out[["N"]] <- n <- nrow(Xy)

  out[["seed"]] <- if(is.null(seed)) sample(10^9, 1) else seed
  set.seed(out$seed)
  if(noisy) message("set seed to ", out$seed, ".\n")
  out[["split"]] <- sample(c("train", "test"), n, replace=TRUE,
                           prob = c(pTraining, pValidation))


  Xy <- as.data.frame(Xy)
  factor_features <- c() # stores individual levels, omitting one

  Xy_distincts <- unlist(lapply(Xy, N_distinct))
  Xy_ch <- unlist(lapply(Xy, is.character))
  for(i in 1:ncol(Xy)){
    if(Xy_distincts[i] == 2 || Xy_ch[i]){
      Xy[,i] <- as.factor(Xy[,i])
    }
  }
  x_factors <- unlist(lapply(Xy[,-ncol(Xy)], is.factor))
  for(i in which(x_factors)){

    if(Xy_distincts[i] > 2){
      tmp <- paste(colnames(Xy)[i], "==", paste0("\'", levels(Xy[,i])[-1], "\'"))
      tmp <- paste0("(", tmp, ")")
      factor_features <- c(factor_features, tmp)
    }else{
      factor_features <- c(factor_features, colnames(Xy)[i])
    }
  }

  # below code can probably
  continuous_features <- cf <- colnames(Xy)[-ncol(Xy)][!x_factors]
  P_continuous <- length(continuous_features)

  out[["continuous_features"]] <- continuous_features
  out[["factors"]] <- colnames(Xy)[-ncol(Xy)][unlist(lapply(Xy[-ncol(Xy)], is.factor))]
  out[["y"]] <- colnames(Xy)[ncol(Xy)]
  out[["P_continuous"]] <- P_continuous
  out[["P_factor"]] <- P_factor <- sum(x_factors)
  P <- P_continuous + P_factor # P does not reflect intercept, interactions, or poly

  if(is.null(outcome)){
    out[["outcome"]] <- if(is.factor(Xy[,ncol(Xy)])) if(N_distinct(Xy[,ncol(Xy)]) > 2) "multinomial" else "binary" else "continuous"
  }

  out[["linear_estimation"]] <- if(out$outcome == "continuous") TRUE else linear_estimation
  if(out$linear_estimation)
    out[["XtX_inv_accepted"]] <- NULL

  out[["y_scale"]] <- if(standardize && out$outcome == "continuous") sd(Xy[out$split == "train", ncol(Xy)]) else 1

  out[["train_scales"]] <- list()
  if(standardize){
    tmp <- which(unlist(lapply(Xy, is_continuous)))
    for(i in tmp){
      out[["train_scales"]][[colnames(Xy)[i]]][["mean"]] <- mean(Xy[,i], na.rm=TRUE)
      out[["train_scales"]][[colnames(Xy)[i]]][["sd"]] <- sd(Xy[,i], na.rm=TRUE)
    }
    Xy[,tmp] <- scale(Xy[,tmp])
  }

  for(i in 2:max_poly_degree)
    continuous_features <- c(continuous_features, paste("pow(", cf, ",", i, ")"))

  features <- f <- c(continuous_features, factor_features)
  if(max_interaction_degree > 1){
    for(i in 2:max_interaction_degree){
      features <- c(f, apply(combn(f, i), 2, paste, collapse = " * "))
    }
  }


  models <- data.frame(features = features,
                       test_adjR2 = NA,
                       stringsAsFactors = FALSE)
  models[["estimated"]]<- FALSE # will be updated based on whether successful ...

  if(out$outcome == "multinomial"){
    out[["best_test_adj_accuracy"]] <- 0
    models[["test_adj_accuracy"]] <- 0
  }else{
    out[["best_test_adjR2"]] <- 0
    # computes pseudo, most important for logit. for ols, multiple definitions of R^2 equivalent including squared correlation.
    # applies the adjustment found in adj R^2 even though out of sample ...
  }

  out[["improvement"]] <- 0  # 0 not meaningful; just initializing ...

  if(out$outcome == "continuous"){
    models[["MAPE"]] <- NA
  }else{
    models[["test_accuracy"]] <- NA
    if(!out$linear_estimation){
      models[["AIC"]] <- NA
      models[["BIC"]] <- NA

    }
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
  out[["threshold_include"]] <- threshold_include
  out[["threshold_estimate"]] <- 0.001
  out[["max_fails"]] <- max_fails
  out[["min_models"]] <- min_models
  out[["file_name"]] <- file_name
  out[["store_fit"]] <- store_fit
  out[["max_block"]] <- 250
  out[["noisy"]] <- noisy

  if(noisy) summary(out, results_overview=FALSE)

  m <- 1            # counts which model

  if(out$outcome == "continuous"){

    y_train <- Xy[out$split == "train", ncol(Xy)]

  }else{ # classification setup

    out[["training_labels"]] <- as.character(unique(Xy[out$split == "train", ncol(Xy)]))
    out[["labels"]] <- Xy[,ncol(Xy)]
    out[["y_test_labels"]] <- out$labels[out$split == "test"]

    if(out$outcome == "multinomial"){

      tallies <- table(Xy[out$split == "train", ncol(Xy)])
      if(noisy && min(tallies) < 10){
        warning("Training the model with rarely observed labels is not recommended.")
        tables(tallies)
      }

      modal_outcome <- names(tallies)[which.max(tallies)]
      out[["reference_category"]] <- modal_outcome
      Xy[ , ncol(Xy)] <- relevel(Xy[,ncol(Xy)], out$reference_category)

      if(out$noisy) message("Multinomial models will be fit with '",
                          modal_outcome,
                          "' (the sample mode of the training data) as the reference category.\n\n")

    }else{

      out[["y_train_mean"]] <- mean(as.numeric(Xy[out$split == "train", ncol(Xy)]) - 1)

    }

    if(out$linear_estimation){

      ln_odds <- log_odds(Xy[,ncol(Xy)], split = (out$split == "train"), noisy = noisy)

      if(out$outcome == "binary"){

        Xy[ , ncol(Xy)] <- ln_odds
        y_train <- ln_odds[out$split == "train"]

      }else{
        y_train <- ln_odds[out$split == "train", ]
      }
    }else{
      y_train <- Xy[out$split == "train", ncol(Xy)]
      y_test <- Xy[out$split == "test", ncol(Xy)]
    }
  } # end classification setup

 if(noisy) message("beginning Forward Stepwise Regression...")

 while((m <= nrow(out$models)) &&
        ((out$improvement > out$threshold_estimate) || m <= out$min_models) &&
        out$unable_to_estimate < out$max_fails){

    out[[mod(m)]] <- list()

    if(sum(out$models$accepted)){

      if(grepl("\\*", out$models$features[m])){

        tmp <- unlist(strsplit(out$models$features[m], "\\*"))
        if(mean(tmp %in% out$models$features[out$models$accepted]) == 1){
          out$models$formula[m] <- paste(out[["best_formula"]], "+", out$models$features[m])
        }else{
          if(noisy) message("Skipping ", out$models$features[m], "\n")
        }
      }else{
        out$models$formula[m] <- paste(out[["best_formula"]], "+", out$models$features[m])
      }
    }else{
      out$models$formula[m] <- paste(out$y_name, "~", out$models$features[m])
    }

    if(!is.na(out$models$formula[m])){

      if(out$outcome == "continuous"){ # fit, etc.

        system.time(out <- ols(out, Xy, m, y = y_train))

      }else{ # begin classification

        if(out$outcome == "multinomial"){

          if(out$linear_estimation){

            system.time(out <- ols(out, Xy, m, y = y_train))

          }else{

            if(noisy) cat("\n\n")
            system.time(out[[mod(m)]][["fit"]] <- multinom(as.formula(out$models$formula[m]),
                                                           Xy[out$split == "train", ], trace = noisy))
            out[[mod(m)]][["coeffs"]] <- t(as.matrix(coefficients(out[[mod(m)]][["fit"]])))
            out <- post_estimation(out, Xy, m)

          }
        }else{ # start binary

          if(out$linear_estimation){
            system.time(out <- ols(out, Xy, m, y = y_train))
          }else{
            system.time(out[[mod(m)]][["fit"]] <- glm(as.formula(out$models$formula[m]),
                                                      Xy[out$split == "train",],
                                                      family = binomial(link = "logit")))
            out[[mod(m)]][["coeffs"]] <- out[[mod(m)]][["fit"]][["coefficients"]]
            out <- post_estimation(out, Xy, m, y_test)
          }
        } # end logit
      } # end fit, etc.

      out$models$estimated[m] <- complete(out[[mod(m)]][["coeffs"]])

      if(out$noisy)
        summary(out, estimation_overview = FALSE, results_overview = FALSE, model_number = m)

      # saving or removing ...
      if(out$store_fit == "none")
        out[[mod(m)]][["fit"]] <- NULL

      if(out$store_fit == "accepted" && !out$models$accepted[m])
        out[[mod(m)]][["fit"]] <- NULL

      if(!is.null(out$file_name)){
        if(out$noisy) message("\n\nsaving (updated) results as ", out$file_name, "\n\n")
        save(out, file=out$file_name)
      }

    }
    m <- m + 1

  } # end WHILE loop

  out$XtX_inv_accepted <- NULL

  if(out$noisy) summary(out, estimation_overview=FALSE)

  return(out)

}


