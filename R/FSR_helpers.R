# if !recursive, divides into blocks based on n and max_block_size
# if recursive, calls block_solve(), rather than solve(), until n/2 < max_block_size
# note: matrix inversion and several matrix multiplications must be performed on largest blocks!
# assumes matrices are dense; otherwise, use sparse options...
# max_block_size chosen by trial-and-error on 2017 MacBook Pro i5 with 16 gigs of RAM
# (too small == too much subsetting, too big == matrix calculations too taxing)
#  S, crossprod(X), will be crossprod(X) only at outer call
# Either S or X should be provided, but not both
# S = | A B |
#     | C D |
# for full expressions used below: https://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion
# returns NULL if inversion fails either due to collinearity or memory exhaustion

block_solve  <- function(S = NULL, X = NULL, max_block_size = 250, A_inv = NULL, recursive=TRUE, noisy=TRUE){

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

  invert <- if(recursive && (k > max_block_size)) block_solve else solvable

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
