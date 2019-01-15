N_distinct <- function(x) if(ncol(as.matrix(x)) == 1) length(unique(x)) else unlist(lapply(x, N_distinct))
#is_continuous <- function(x) if(is.numeric(x)) N_distict(x) > 2 else FALSE
is_continuous <- function(x) unlist(lapply(x, is.numeric)) & N_distinct(x) > 2
mod <- function(m) paste0("model", m)
complete <- function(x) !is.null(x) && sum(is.na(x)) == 0
match_arg <- function(arg, choices){if(is.null(arg)) arg else match.arg(arg, choices)}

# ensure unique column names with pow
# degree may be vector
pow <- function(X, degree){

  if(length(degree) == 1){
    X <- as.matrix(X)^degree
    colnames(X) <- paste0(colnames(X), "_deg_", degree)
  }else{
    out <- list()
    for(i in 1:length(degree))
      out[[i]] <- as.matrix(pow(X, degree[i]))
    X <- do.call(cbind, out)
  }
  return(X)

}
`%^%` <- function(X, degree) pow(X, degree)

model_matrix <- function(f, d, intercept, noisy=TRUE){

  tried <- try(model.matrix(as.formula(f), d, na.action = "na.omit"), silent=TRUE)

  if(inherits(tried, "try-error")){
    if(noisy) cat("model.matrix() reported the following error:\n", tried, "\n\n")
    return(NULL)
  } else {

    attr(tried, which = "formula") <- f # works in local but not global environment, fix

    if(intercept) return(tried) else return(tried[,-1])
  }

}

get_interactions <- function(features, maxInteractDeg, may_not_repeat = NULL){

  interactions <- list()

  for(i in 2:maxInteractDeg){

    combos <- combn(features, i) # i x choose(n, i) matrix
    combos <- combos[ , which_include(combos, may_not_repeat)]

    interactions[[i]] <- apply(combos, 2, paste, collapse = " * ")

  }

  return(unlist(interactions))

}

which_include <- function(combos, may_not_repeat){
# prevents multiplication of mutually exclusive categorical variables' levels
# suppose you have a factor variable, party with levels D, R, I
# at this point, factor features are strings formatted
# (party == 'D') and (party == 'R')
# but identical((party == 'D') * (party == 'R'), rep(0, N)) == TRUE
# this function uses grepl() to prevent such 0 columns from entering
# the formula subsequently...
#
# also, different monomials of the same variable should not interact
# raising the polynomial degree beyond user specification

  keepers <- 1:ncol(combos)

  if(length(may_not_repeat) == 0){

    return(keepers)

  }else{

    to_drop <- list()

    for(i in 1:length(may_not_repeat)){
      to_drop[[i]] <- which(colSums(apply(combos, 2, grepl, pattern = may_not_repeat[i])) > 1)
    }
    to_drop <- unique(unlist(to_drop))

    if(length(to_drop)) return(keepers[-to_drop]) else return(keepers)
  }
}

# depracated
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

classify <- function(probs, as_factor=TRUE, labels=NULL, cutoff = NULL){ # not meant for binary labels...

  if(ncol(as.matrix(probs)) == 1){

    if(is.null(labels))
      labels <- c("label1", "label2")

    if(is.null(cutoff))
      cutoff <- 0.5

    classified <- labels[(probs > cutoff) + 1]

  }else{

    if(!is.null(labels))
      colnames(probs) <- labels
    if(is.null(colnames(probs)))
      colnames(probs) <- paste0("label", 1:ncol(probs))
    classified <- colnames(probs)[apply(probs, 1, which.max)]

  }

  if(as_factor)
    classified <- as.factor(classified)

  return(classified)

}


log_odds <- function(x, split = NULL, noisy = TRUE){

  if(N_distinct(x) == 2){

    if(is.factor(x))
      x <- as.numeric(x) - 1
    p <- mean(if(is.null(split)) x else x[split], na.rm=TRUE)
    y <- ifelse(x == 1, log(p/(1 - p)), log((1 - p)/p))

  }else{

    if(!is.factor(x))
      x <- as.factor(x)

    if(is.null(split)){

      p_reference <- mean(x == levels(x)[1])
      y <- matrix(nrow = length(x), ncol = (length(levels(x)) - 1))
      colnames(y) <- levels(x)[-1]
      for(i in 1:ncol(y)){
        p_interest <- mean(x == levels(x)[i + 1])
        y[ , i] <- ifelse(x == levels(x)[i + 1],
                          log(p_interest/p_reference),
                          log(p_reference/p_interest))
      }

    }else{ # put whole sample on training scale, so N rows, not N_train

      x_train <- x[split]
      p_reference <- mean(x_train == levels(x)[1])
      y <- matrix(nrow = length(x),
                  ncol = (length(levels(x_train)) - 1))
      colnames(y) <- levels(x_train)[-1]
      for(i in 1:ncol(y)){
        p_interest <- mean(x_train == levels(x_train)[i + 1])
        y[ , i] <- ifelse(x == levels(x_train)[i + 1],
                          log(p_interest/p_reference),
                          log(p_reference/p_interest))
      }
    }
  }

  if(noisy && sum(is.na(y)))
    warning("NAs encountered by log_odds")

  return(y)

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



ols <- function(object, Xy, m, train = TRUE, y = NULL, y_test = NULL){

  X <- if(train){
          model_matrix(formula(object$models$formula[m]),
                        Xy[object$split == "train", ],
                        noisy = object$noisy, intercept=TRUE)
       }else{
          model_matrix(formula(object$models$formula[m]),
                        Xy, noisy = object$noisy, intercept=TRUE)
       }

  if(exists("X")){

    if(is.null(y))
      y <- if(train) Xy[object$split == "train", ncol(Xy)] else Xy[, ncol(Xy)]

    if(ncol(X) >= length(y) && object$noisy){

      message("There are too few training observations to estimate further models (model == ",
              m, "). Exiting.")
      object$unable_to_estimate <- object$max_fails

    }else{

      XtX_inv <- block_solve(X = X, max_block = object$max_block,
                             A_inv = object$XtX_inv_accepted)
                                    # initialized to NULL,
                                    # which block_solve interprets as 'start from scratch'

      if(!is.null(XtX_inv)){

        object[[mod(m)]][["coeffs"]] <- tcrossprod(XtX_inv, X) %*% y

        if(complete(object[[mod(m)]][["coeffs"]])){

          object$models$estimated[m] <- TRUE

          object <- post_estimation(object, Xy, m, y_test)
          if(object$models$accepted[m])
            object$XtX_inv_accepted <- XtX_inv

          remove(XtX_inv)

        }
      }
    }
  }
  if(!object$models$estimated[m]){
    warning("Unable to estimate model", m, "\n\n")
    object$unable_to_estimate <- object$unable_to_estimate + 1
  }
  if(object$noisy) cat("\n")
  return(object)
}

post_estimation <- function(object, Xy, m, y_test = NULL){

    P <- if(object$outcome == "multinomial")
            nrow(object[[mod(m)]][["coeffs"]]) else length(object[[mod(m)]][["coeffs"]])

    object$models$P[m] <- object[[mod(m)]][["p"]]  <- P

    if(is.null(y_test))
      y_test <- Xy[object$split == "test", ncol(Xy)]

    if(object$outcome == "continuous"){

      object[[mod(m)]][["y_hat"]] <- predict(object, Xy[object$split=="test", ], m, standardize = FALSE)
      MAPE <- object$y_scale * mean(abs(object[[mod(m)]][["y_hat"]] - y_test))
      object$models$MAPE[m] <- object[[mod(m)]][["MAPE"]] <- MAPE

    }else{

      pred <- predict(object, Xy[object$split=="test", ], m, standardize = FALSE)

      object[[mod(m)]][["y_hat"]] <- pred$probs
      object[[mod(m)]][["classified"]] <- pred$classified

      object$models$test_accuracy[m] <- mean(as.character(pred$classified) == object$y_test_labels)

      if(!object$linear_estimation){

        object$models$AIC[m] <- if(object$outcome == "binary")
                                    object[[mod(m)]][["fit"]][["aic"]] else
                                      object[[mod(m)]][["fit"]][["AIC"]]

        object$models$BIC[m] <- object$models$AIC[m] - 2*P + log(object$N_train)*P

      }
    }

    if(object$outcome != "multinomial"){

      R2 <- cor(object[[mod(m)]][["y_hat"]], as.numeric(y_test))^2

      adjR2 <- (object$N_train - P - 1)/(object$N_train - 1)*R2

      object$models$test_adjR2[m] <- object[[mod(m)]][["adj_R2"]] <- adjR2

      improvement <- adjR2 - object$best_test_adjR2

    }else{

      adj_accuracy <- (object$N_train - P)/(object$N_train - 1)*object$models$test_accuracy[m]

      object$models$test_adj_accuracy[m] <- adj_accuracy

      improvement <- adj_accuracy - object$best_test_adj_accuracy

    }

    object[["improvement"]] <- improvement

  if(object$improvement > object$threshold_include){

      object[["best_formula"]] <- object$models$formula[m]
      object[["best_coeffs"]] <- object[[mod(m)]][["coeffs"]]

      if(object$outcome == "multinomial"){
        object[["best_test_adj_accuracy"]] <- adj_accuracy
      }else{
        object[["best_test_adjR2"]] <- adjR2
      }

      object$models$accepted[m] <- TRUE

      if(object$outcome == "continuous")
        object[["best_MAPE"]] <- MAPE
  }
  return(object)
}


