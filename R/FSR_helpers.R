FSR_estimate <- function(z, Xy){

  y_train <- Xy[z$split == "train", ncol(Xy)]
  y_test <- Xy[z$split == "test", ncol(Xy)]

  m <- 1            # counts which model
  improvement <- 0  # 0 not meaningful; just initializing ...

  while((m <= nrow(z$models)) &&
        ((improvement > z$threshold_estimate) || m <= z$min_models) &&
        z$unable_to_estimate < z$max_fails){

    if(z$noisy) cat("\n\n\n\nModel:", m, "\n\n")
    z[[mod(m)]] <- list()
    z$models$formula[m] <- if(sum(z$models$accepted) == 0) paste(z$y_name, "~", z$models$features[m]) else paste(z[["best_formula"]], "+", z$models$features[m])

    if(z$model_type == "lm"){

      X_train <- model_matrix(formula(z$models$formula[m]), Xy[z$split == "train",], noisy=z$noisy)

      if(!exists("X_train")){
        if(z$noisy) message("Unable to construct model.matrix for model ",  m, ". Skipping.")
        z$unable_to_estimate <- z$unable_to_estimate + 1
      }else{

        z$models$p[m] <- z[[mod(m)]][["p"]] <- ncol(X_train)

        if(z[[mod(m)]][["p"]] >= z$N_train){
          if(z$noisy) message("There are too few training observations to estimate model ",  m, ". Skipping.")
          z$unable_to_estimate <- z$unable_to_estimate + 1
        }else{

          if(sum(z$models$accepted) == 0){ # passing X takes crossproduct first; starting with second iteration,
            XtX_inv <- block_solve(X = X_train, max_block = z$max_block)
          }else{  # XtX_inv of the last accepted model is taken as the inverse of the first block
            XtX_inv <- block_solve(X = X_train, A_inv = XtX_inv_accepted, max_block = z$max_block)
          }
          if(!is.null(XtX_inv)){

            z[[mod(m)]][["coeffs"]] <- tcrossprod(XtX_inv, X_train) %*% y_train

            if(m == length(z$models$features)) remove(XtX_inv)

            if(!complete(z[[mod(m)]][["coeffs"]])){

              if(z$noisy) cat("\nfailed to estimate model", m, "skipping...\n")
              z$unable_to_estimate <- z$unable_to_estimate + 1

            }else{

              z[[mod(m)]][["y_hat"]] <- model_matrix(formula(z$models$formula[m]),
                                                     Xy[z$split == "test", ]) %*% z[[mod(m)]][["coeffs"]]

              R2 <- cor(z[[mod(m)]][["y_hat"]], y_test, method=z$cor_type)^2
              adjR2 <- (z$N_train - ncol(X_train) - 1)/(z$N_train - 1)*R2 # odd to have penalty mashed up this way...
              z$models$adjR2[m] <- z[[mod(m)]][[paste0("adj_R2_", z$cor_type)]] <- adjR2

              MAPE <- z$y_scale * mean(abs(z[[mod(m)]][["y_hat"]] - y_test))
              z$models$MAPE[m] <- z[[mod(m)]][["MAPE"]] <- MAPE
              improvement <- if(sum(z$models$accepted) == 0) adjR2 else (adjR2 - z$best_adjR2)

              if(improvement > z$threshold_include){
                z$best_formula <- z$models$formula[m]
                z[["best_coeffs"]] <- z[[mod(m)]][["coeffs"]]
                z[["best_adjR2"]] <- adjR2
                z[["best_MAPE"]] <- MAPE
                z$models$accepted[m] <- TRUE
                XtX_inv_accepted <- XtX_inv
              }

              if(z$noisy) summary(z, estimation_overview=FALSE, results_overview=FALSE, model_number = m)
            }
          }else{
            if(z$noisy) cat("\nunable to estimate Model", m, "likely due to (near) singularity.\n")
            unable_to_estimate <- unable_to_estimate + 1
          }
        }
      }# end linear model

    }else{ # begin classification

      if(m == 1)
        z[["training_labels"]] <- as.character(unique(Xy[z$split == "train", ncol(Xy)])) # levels(Xy[[z$split]])

      if(z$model_type == "multinomial"){

        if(m == 1){

          tallies <- table(Xy[z$split == "train", ncol(Xy)])
          modal_outcome <- names(tallies)[which.max(tallies)]
          z[["reference_category"]] <- modal_outcome
          Xy[,ncol(Xy)] <- relevel(Xy[,ncol(Xy)], z$reference_category)

          if(z$noisy) message("Multinomial models will be fit with '",
                            modal_outcome,
                            "' (the sample mode of the training data) as the reference category.\n\n")

          z[["y_train_means"]] <- tallies/z$N_train # proportions... may not use...

        }

        z[[mod(m)]][["fit"]] <- multinom(as.formula(z$models$formula[m]),
                                         Xy[z$split == "train", ])
        z[[mod(m)]][["fit"]][["aic"]] <- z[[mod(m)]][["fit"]][["AIC"]]


        z[[mod(m)]][["coeffs"]] <- coefficients(z[[mod(m)]][["fit"]])
        z$models$P[m] <- z[[mod(m)]][["p"]] <- ncol(z[[mod(m)]][["coeffs"]])

        # defaults to mean-bias reducing adjusted scores, see ?brglmFit for 'type' options including ML
        # z[[mod(m)]][["fit"]] <- brmultinom(as.formula(z$models$formula[m]),
        #                                   Xy[z$split == "train", ],
        #                                   ref = z$reference_category)
        # for training accuracy, could precede as below
        # fitted values are predicted probabilities for reference category, followed by others, as vector...
        # z[[mod(m)]][["probs"]] <- matrix(z[[mod(m)]][["fit"]][["fitted.values"]], ncol=length(z$training_labels))

        # classified <- classify(z[[mod(m)]][["probs"]],
        #                       labels = c(z$reference_category,
        #                                  z$training_labels[-which(z$training_labels == z$reference_category)]))


        predictions <- predict(z, Xy[z$split == "test", ], model_to_use = m, standardize = FALSE)
        z[[mod(m)]][["pred_probs"]] <- predictions$probs
        z[[mod(m)]][["classified"]] <- predictions$classified

        z$models$test_accuracy[m] <- mean(y_test == as.character(z[[mod(m)]][["classified"]]))
        z$models$test_adj_accuracy[m] <- adj_accuracy <- (z$N_train - z[[mod(m)]][["p"]])/(z$N_train - 1)*z$models$test_accuracy[m]
        improvement <- if(sum(z$models$accepted)) (adj_accuracy - z$best_adj_accuracy) else adj_accuracy

        if(improvement > z$threshold_include){

          z[["best_adj_accuracy"]] <- adj_accuracy

        }


      }else{ # start logistic regression code

        if(m == 1){
          z[["y_train_mean"]] <- mean(as.numeric(Xy[z$split == "train",ncol(Xy)]) - 1)
        }

        z[[mod(m)]][["fit"]] <- glm(as.formula(z$models$formula[m]), Xy[z$split == "train",],
                                    family = binomial(link = "logit"))
        z[[mod(m)]][["coeffs"]] <- beta_hat <- z[[mod(m)]][["fit"]][["coefficients"]]

        z$models$p[m] <- z[[mod(m)]][["p"]] <- length(beta_hat)

        z[[mod(m)]][["y_hat"]] <- predict(z[[mod(m)]][["fit"]],
                                          as.data.frame(model_matrix(as.formula(z$models$formula[m]),
                                                                     Xy[z$split == "test",])))

        pseudo_R2 <- cor(z[[mod(m)]][["y_hat"]], as.numeric(y_test))^2
        z$models$adjR2[m] <- adjR2 <- (z$N_train - z[[mod(m)]][["p"]])/(z$N_train - 1)*pseudo_R2

        z[[mod(m)]][["classified"]] <- factor(z[["training_labels"]][(z[[mod(m)]][["y_hat"]] > z[["y_train_mean"]]) + 1],
                                              levels = levels(y_train))
        z$models$test_accuracy[m] <- mean(z[[mod(m)]][["classified"]] == y_test)

        improvement <- if(sum(z$models$accepted)) (adjR2 - z$best_adjR2) else adjR2

        if(improvement > z$threshold_include){

          z[["best_adjR2"]] <- adjR2

        }

      } # end logit

      if(improvement > z$threshold_include){

        z$best_formula <- z$models$formula[m]
        z[["best_coeffs"]] <- z[[mod(m)]][["coeffs"]]
        z$models$accepted[m] <- TRUE
        z[["best_test_accuracy"]] <- z$models$test_accuracy[m]

      }

      z$models$AIC[m] <- z[[mod(m)]][["fit"]][["aic"]]
      z$models$BIC[m] <- z$models$AIC[m] - 2*z[[mod(m)]][["p"]]  + log(z$N_train)*z[[mod(m)]][["p"]]
      if(z$noisy) summary(z, estimation_overview=FALSE, results_overview=FALSE, model_number = m)

      }

    z$models$estimated[m] <- complete(z[[mod(m)]][["coeffs"]])
    if(z$store_fit == "none")
      z[[mod(m)]][["fit"]] <- NULL
    if(z$store_fit == "accepted" && !z$models$accepted[m])
      z[[mod(m)]][["fit"]] <- NULL

    m <- m + 1

    if(!is.null(z$file_name)){
      if(z$noisy) message("\n\nsaving (updated) results as ", z$file_name, "\n\n")
      save(z, file=z$file_name)
    }

  } # end WHILE loop

  if(z$noisy) summary(z, estimation_overview=FALSE)

  if(sum(z$models$estimated) == 0) return(NULL) else return(z)

}


# if !recursive, divides into blocks based on n and max_block
# if recursive, calls block_solve(), rather than solve(), until n/2 < max_block
# note: matrix inversion and several matrix multiplications must be performed on largest blocks!
# assumes matrices are dense; otherwise, use sparse options...
# max_block chosen by trial-and-error on 2017 MacBook Pro i5 with 16 gigs of RAM
# (too small == too much subsetting, too big == matrix calculations too taxing)
#  S, crossprod(X), will be crossprod(X) only at outer call
# Either S or X should be provided, but not both
# S = | A B |
#     | C D |
# for full expressions used below: https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion
# returns NULL if inversion fails either due to collinearity or memory exhaustion

block_solve  <- function(S = NULL, X = NULL, max_block = 250, A_inv = NULL, recursive=TRUE, noisy=TRUE){

  if(is.null(S) == is.null(X))
    stop("Please provide either rectangular matrix as X or a square matrix as S to be inverted by block_solve(). (If X is provided, (X'X)^{-1} is returned but in a more memory efficient manner than providing S = X'X directly).")

  if(!is.null(A_inv) && is.null(X))
    stop("If A_inv is provided, X must be provided to block_solve() too. (Suppose A_inv has p columns; A must be equal to solve(crossprod(X[,1:p])) or, equivalently, block_solve(X=X[,1:p]).")

  solvable <- function(A, noisy=TRUE){

    tried <- try(solve(A), silent = noisy)
    if(noisy) cat(".")
    if(inherits(tried, "try-error")) return(NULL) else return(tried)

  }

  if(is.null(X)){

    stopifnot(nrow(S) == ncol(S))

    symmetric <- isSymmetric(S)
    n <- ncol(S)   # if S is crossprod(X), this is really a p * p matrix
    k <- floor(n/2)

    A <- S[1:k, 1:k]
    B <- S[1:k, (k + 1):n]
    D <- S[(k + 1):n, (k + 1):n]

  }else{

    n <- ncol(X)     # n refers to the resulting crossproduct of S as above
    if(is.null(A_inv)){
      k <- floor(n/2)
      A <- crossprod(X[,1:k])
    }else{
      k <- ncol(A_inv)
    }
    B <- crossprod(X[,1:k], X[,(k+1):n])
    D <- crossprod(X[,(k+1):n])

    symmetric <- TRUE   # refers to S, not A, B, or D (B in general will be rectangular...)

  }

  invert <- if(recursive && (k > max_block)) block_solve else solvable

  if(is.null(A_inv)){
    A_inv <- invert(A, noisy=noisy)
    remove(A)
  }

  if(!is.null(A_inv)){


    if(symmetric){
      # S, crossprod(X), will be symmetric at highest level but not at lower levels
      # want memory savings from that symmetry when it applies
      # by symmetry, B == t(C), so C is never constructed
      if(exists("S")) remove(S)
      C.A_inv <- crossprod(B, A_inv) # really C %*% A_inv since C == t(B)
      schur_inv <- invert(D - C.A_inv %*% B)
      remove(D)

      if(!is.null(schur_inv)){

        S_inv <- matrix(nrow=n, ncol=n)

        S_inv[1:k, 1:k] <- A_inv + A_inv %*% B %*% schur_inv %*% C.A_inv
        remove(B, A_inv)
        S_inv[(k+1):n, 1:k] <- -schur_inv %*% C.A_inv
        S_inv[(k+1):n, (k+1):n] <- schur_inv
        remove(schur_inv, C.A_inv)
        S_inv[1:k, (k+1):n] <- t(S_inv[(k+1):n, 1:k]) # since symmetric matrices have symm inverses
        return(S_inv)

      }else{
        return(NULL)
      }

    }else{

      C.A_inv <- crossprod(B, A_inv) # S[(k+1):n, 1:k] %*% A_inv  # really C %*% A_inv

      if(exists("C.A_inv")){

        if(exists("S")) remove(S)

        schur_inv <- invert(D - C.A_inv %*% B, noisy=noisy)
        remove(D)

        S_inv <- matrix(nrow=n, ncol=n)
        S_inv[1:k, 1:k] <- A_inv + A_inv %*% B %*% schur_inv %*% C.A_inv
        S_inv[(k+1):n, 1:k] <- -schur_inv %*% C.A_inv
        remove(C.A_inv)
        S_inv[(k+1):n, (k+1):n] <- schur_inv
        S_inv[1:k, (k+1):n] <- -A_inv %*% B %*% schur_inv
        remove(B, A_inv, schur_inv)
        return(S_inv)

      }else{
        return(NULL)
      }
    }
  }else{
    return(NULL)
  }

}

N_distinct <- function(x) length(unique(x))
is_continuous <- function(x) if(is.numeric(x)) length(unique(x)) > 2 else FALSE
mod <- function(m) paste0("model", m)

pow <- function(X, degree){

  X <- as.matrix(X)^degree
  colnames(X) <- paste0(colnames(X), "_deg_", degree)
  return(X)    # ensure unique column names

}

model_matrix <- function(f, d, noisy=TRUE){
  tried <- try(model.matrix(f, d), silent=TRUE)
  if(inherits(tried, "try-error")){
    if(noisy) cat("model.matrix() reported the following error:\n", tried, "\n\n")
    return(NULL)
  } else {
    return(tried)
  }
}

isolate_interaction <- function(elements, degree){

  f <- paste(elements, collapse = " * ")
  for(i in 1:degree){
    tmp <- combn(elements, i)
    if(i > 1)
      tmp <- apply(tmp, 2, paste, collapse="*")
    f <- paste(f, "-", paste(tmp, collapse=" - "))
  }
  return(f)

}

classify <- function(probs, as_factor=TRUE, labels=NULL){ # not meant for binary labels...

  if(!is.null(labels))
    colnames(probs) <- labels
  if(is.null(colnames(probs)))
    colnames(probs) <- paste0("label", 1:ncol(probs))
  classified <- colnames(probs)[apply(probs, 1, which.max)]
  if(as_factor)
    classified <- as.factor(classified)
  return(classified)

}

complete <- function(x){
  !is.null(x) && sum(is.na(x)) == 0
}

