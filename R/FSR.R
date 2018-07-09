block_solve  <- function(S = NULL, X = NULL, max_block_size = 250, recursive=TRUE, noisy=TRUE){

  # if !recursive, divides into blocks based on n and max_block_size
  # if recursive, calls block_solve(), rather than solve(), until n/2 < max_block_size

  # note: matrix inversion and several matrix multiplications must be performed on largest blocks!
  # assumes matrices are dense; otherwise, use sparse options...

  # max_block_size chosen by trial-and-error on 2017 MacBook Pro i5 with 16 gigs of RAM
  # (too small == too much subsetting, too big == matrix calculations too taxing)

  #  S, crossprod(X), will be crossprod(X) at outer call

  # Either S or X should be provided, but not both

  # S = | A B |
  #     | C D |
  # for full expressions used below: https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion


  if(is.null(X)){

    stopifnot(nrow(S) == ncol(S))  # called internally, may not be necessary checks

    symmetric <- isSymmetric(S)
    n <- ncol(S)   # if S is crossprod(X), this is really a p * p matrix
    k <- floor(n/2)

    A <- S[1:k, 1:k]
    B <- S[1:k, (k + 1):n]
    D <- S[(k + 1):n, (k + 1):n]

  }else{

    n <- ncol(X)     # n refers to the resulting crossproduct of S as above
    k <- floor(n/2)

    A <- crossprod(X[,1:k])
    B <- crossprod(X[,1:k], X[,(k+1):n])
    D <- crossprod(X[,(k+1):n])

    symmetric <- TRUE   # refers to S, not A, B, or D (B in general will be rectangular...)

  }

  solvable <- function(A, noisy=noisy){

    tried <- try(solve(A), silent = TRUE)
    if(inherits(tried, "try-error")) return(NULL) else return(tried)

  }

  invert <- if(recursive && (k > max_block_size)) block_solve else solvable

  A_inv <- invert(A)
  remove(A)

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

        schur_inv <- invert(D - C.A_inv %*% B)
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
# returns list with slope coefficients, model details, and fit vs. the test data of the final model
#' @export
FSR <- function(Xy,
                max_poly_degree = 4, max_interaction_degree = 2,
                cor_type = "kendall", threshold = 0.001,
                pTraining = 0.7, pValidation = 0.2, pTesting = 0.1, max_block_size = 250,
                noisy = TRUE, seed = NULL,
                model = "lm"){

  N_distinct <- function(x) length(unique(x))
  is_continuous <- function(x) if(is.numeric(x)) length(unique(x)) > 2 else FALSE
  mod <- function(m) paste0("model", m)
  pow <- function(X, degree){
    X <- X^degree
    colnames(X) <- paste0(colnames(X), "_deg_", degree)
    return(X)    # ensure unique column names
  }

  if(!is.matrix(Xy) && !is.data.frame(Xy))
    stop("Xy must be a matrix or data.frame. Either way, y must be the final column.")
  if(min(pTraining, pTesting, pValidation) < 0 || max(pTraining, pTesting, pValidation) > 1)
    stop("pTraining, pTesting, and pValidation should all be between 0 and 1 and sum to 1.")
  stopifnot(is.numeric(threshold))

  n <- nrow(Xy)

  Xy <- as.data.frame(Xy)
  N_factor_columns <- 0
  for(i in 1:ncol(Xy)){
    k <- N_distinct(Xy[,i])
    if(k == 2 || is.character(Xy[,i])){
        Xy[,i] <- as.factor(Xy[,i])    # switch this to as.character() in case of nuissance
        N_factor_columns <- N_factor_columns + k - 1
    }
  }
  continuous_features <- colnames(Xy)[-ncol(Xy)][unlist(lapply(Xy[-ncol(Xy)], is_continuous))]
  P_features <- length(continuous_features) + N_factor_columns # columns without intercept...

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
      warning("Logistic regression not implemented, treating as continuous (for now...)...")
    }
  }

  out <- list()
  out[["seed"]] <- if(is.null(seed)) as.integer(runif(1, 0, 10000000)) else seed
  set.seed(out$seed)
  if(noisy) message("set seed to ", out$seed, ".\n")

  out[["split"]] <- sample(c("train", "validate", "test"), n, replace=TRUE,
                  prob = c(pTraining, pValidation, pTesting))
  y_train <- Xy[out$split == "train", ncol(Xy)]
  y_validate <- Xy[out$split == "validate", ncol(Xy)]
  y_test <- Xy[out$split == "test", ncol(Xy)]

  increment <- "neither"
  for(i in 1:min(max_poly_degree, max_interaction_degree))
    increment <- c(increment, "poly", "interaction")
  increment <- c(increment, rep("poly", max_poly_degree - min(max_poly_degree, max_interaction_degree)))
  increment <- c(increment, rep("interaction", max_interaction_degree - min(max_poly_degree, max_interaction_degree)))
  out[["estimated"]] <- rep(FALSE, length(increment)) # will be updated based on whether successful ...

  m <- 1            # counts which model
  improvement <- threshold + 1  # not meaningful; just initializing ...
  poly_degree <- 1
  interaction_degree <- 0
  out[[mod(m)]][["formula"]] <- formula(paste(colnames(Xy)[ncol(Xy)], "~ ."))
  N_train <- sum(out$split == "train")

  while((m <= length(increment)) && (improvement > threshold) && (P_features < N_train)){

    if(increment[m] == "poly"){

      poly_degree <- poly_degree + 1
      Xy <- cbind(pow(Xy[ , match(continuous_features, colnames(Xy))], poly_degree), Xy)
      out[[mod(m)]][["formula"]] <- out[[mod(m - 1)]][["formula"]]

    }

    if(increment[m] == "interaction"){

      f <- paste(colnames(Xy)[ncol(Xy)], "~ .")
      f <- paste(c(f, rep("*.", sum(increment[1:m] == "interaction"))), collapse="")
      out[[mod(m)]][["formula"]] <- formula(f)

    }

    X_train <- model.matrix(out[[mod(m)]][["formula"]], Xy[out$split == "train", ])
    out[[mod(m)]][["p"]] <- ncol(X_train)

    if(out[[mod(m)]][["p"]] >= nrow(X_train)){

      if(noisy) message("There are too few training observations to estimate model ",  m, ". Skipping.")

    }else{

      XtX_inv <- block_solve(X = X_train,  # passing X takes crossproduct first
                             max_block_size)

      if(!is.null(XtX_inv)){

        out[[mod(m)]][["est"]] <- tcrossprod(XtX_inv, X_train) %*% y_train
        remove(XtX_inv)

        out[[mod(m)]][["y_hat"]] <- model.matrix(out[[mod(m)]][["formula"]],
                                                 Xy[out$split == "validate", ]) %*% out[[mod(m)]][["est"]]

        out[[mod(m)]][[paste0("R2_", cor_type)]] <- cor(out[[mod(m)]][["y_hat"]], y_validate, method=cor_type)^2

        improvement <- if(m == 1) out[[mod(m)]][[paste0("R2_", cor_type)]] else out[[mod(m)]][[paste0("R2_", cor_type)]] - out[[mod(m - 1)]][[paste0("R2_", cor_type)]]

        out$estimated[m] <- TRUE
        if(m == 1 && noisy) cat("\tpseudo R^2 (", cor_type, ")\n", sep="")
        if(noisy) cat("Model ", m, ": ", out[[mod(m)]][[paste0("R2_", cor_type)]], "\n", sep="")

      }else{
        if(noisy) cat("unable to estimate Model", m, "\n")
      }
    }

    m <- m + 1
    P_features <- P_features + (increment[m] == "poly") * length(continuous_features) + (increment[m] == "interaction") * choose(P_features, interaction_degree + 1)

  } # end WHILE loop

  if(sum(out$estimated) == 0){

    if(noisy) message("No models could be estimated, likely due to (near) singularity; returning NULL. Check for highly correlated features or factors with rarely observed levels. (Increasing pTraining may help.)\n")
    return(NULL)

  }else{

    # X_test used only once; if performance is much worse than on X_validate, sign of overfitting

    if(sum(out$estimated) == 1){

      selected <- which(out$estimated)

    }else{

      last_fit <- max(which(out$estimated))
      second_to_last <- max(which(out$estimated[-last_fit]))
      selected <- if(out[[mod(last_fit)]][[paste0("R2_", cor_type)]] > out[[mod(second_to_last)]][[paste0("R2_", cor_type)]]) last_fit else second_to_last
      if(sum(increment[1:selected] == "poly") < sum(increment[1:m] == "poly"))
        Xy <- Xy[, -match(paste0(continuous_features, "_deg_", sum(increment[1:m] == "poly") + 1), colnames(Xy))]

    }

    out[["test"]][["y_hat"]] <- model.matrix(out[[mod(selected)]][["formula"]],
                                             Xy[out$split == "test", ]) %*% out[[mod(selected)]][["est"]]
    out[["test"]][[paste0("R2_", cor_type)]] <- cor(out[["test"]][["y_hat"]], y_test, method=cor_type)^2
    out$MAPE <- mean(abs(out$test$y_hat - y_test))

    if(noisy) cat("\n\nModel ", selected,
                  " (polynomial degree ", sum(increment[1:selected] == "poly") + 1,
                  ", interactions of order ", sum(increment[1:selected] == "interactions"),
                  ") was chosen since it had the highest Pseudo R^2 (", cor_type,
                  ") on the validation data (",
                  out[[mod(selected)]][[paste0("R2_", cor_type)]],
                  "). The model has Mean Absolute Predicted Error on the test data of ",
                  out$MAPE, " (Pseudo R^2 ",  cor_type,
                  " on the test data: ", out[["test"]][[paste0("R2_", cor_type)]], ").", sep="")

    return(out)
  }
}
