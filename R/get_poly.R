get_poly <- function(xdata = NULL, deg, maxInteractDeg = deg, Xy = NULL, warn_level=2){

 # if(is.null(xdata) && is.null(Xy))
 #   stop("please provide either xdata or Xy (if provided, Xy must have the dependent variable in the final column).")

  # helper functions
  N_distinct <- function(x) length(unique(x))
  is_continuous <- function(x) if(is.numeric(x)) length(unique(x)) > 2 else FALSE
  mod <- function(m) paste0("model", m)
  pow <- function(X, degree){
    X <- X^degree
    colnames(X) <- paste0(colnames(X), "_deg_", degree)
    return(X)    # ensure unique column names, necessary for model.matrix()
  }
  model_matrix <- function(f, d){
    tried <- try(model.matrix(f, d), silent=TRUE)
    if(inherits(tried, "try-error")){
      cat("model.matrix() reported the following error:\n", tried, "\n\n")
      return(NULL)
    } else {
      return(tried)
    }
  }

  if(!is.null(xdata)){

    if(!is.matrix(xdata) && !is.data.frame(xdata))
      stop("xdata must be a matrix or data.frame.")

    N <- nrow(xdata)
    xdata <- as.data.frame(xdata)

    N_factor_columns <- 0
    for(i in 1:ncol(xdata)){
      k <- N_distinct(xdata[,i])
      if(is.character(xdata[,i]) || k == 2){
        xdata[,i] <- as.factor(xdata[,i])
        N_factor_columns <- N_factor_columns + k - 1
      }
    }

    f <- "~ ."

    if(deg > 1){
      which_continuous <- which(unlist(lapply(xdata, is_continuous)))
      for(d in 2:deg)
        xdata <- cbind(xdata, pow(xdata[ , which_continuous], d))
    }

    P <- sum(choose(N_factor_columns + deg*sum(which_continuous), 1:maxInteractDeg)) + 1

    if(P > N){
      if(warn_level == 2){
        message("The model matrix will have ", N, " observations and ", P, " columns.\n")
        stop("\n\ntoo few observations to estimate. set poly_reg(..., warn_level=1) to warn without stopping or poly_reg(..., warn_level=0) to bypass (e.g., to create X_test).")
      }
      if(warn_level == 1){
        message("The model matrix will have ", N, " observations and ", P, " columns.\n")
        warning("\n\ntoo few observations to estimate. set poly_reg(..., warn_level=2) to stop on this error or poly_reg(..., warn_level=0) to bypass.")
      }
    }



    if(maxInteractDeg > 0)
      f <- paste(c(f, rep("*.", maxInteractDeg)), collapse="")

    return(model_matrix(formula(f), xdata))

  }else{

    # this will only return a different object than the above if y contains missing data

    if(!is.matrix(Xy) && !is.data.frame(Xy))
      warning("Xy must be a matrix or data.frame. Either way, y must be the final column.")

    N <- nrow(Xy)
    Xy <- as.data.frame(Xy)

    N_factor_columns <- 0
    for(i in 1:ncol(Xy)){
      k <- N_distinct(Xy[,i])
      if(is.character(Xy[,i]) || k == 2){
        Xy[,i] <- as.factor(Xy[,i])
        N_factor_columns <- N_factor_columns + k - 1
      }
    }

    f <- paste(colnames(Xy)[ncol(Xy)], "~ .")

    if(deg > 1){
      which_continuous <- which(unlist(lapply(Xy[-ncol(Xy)], is_continuous)))
      for(d in 2:deg)
        Xy <- cbind(Xy, pow(Xy[ , which_continuous], d))
    }

    P <- sum(choose(N_factor_columns + deg*sum(which_continuous), 1:maxInteractDeg)) + 1

    if(P > N){
      if(warn_level == 2){
        message("\n\nThe model matrix will have ", N, " observations and ", P, " columns.\n")
        stop("too few observations to estimate. set poly_reg(..., warn_level=1) to warn without stopping or poly_reg(..., warn_level=0) to bypass.")
      }
      if(warn_level == 1){
        message("\n\nThe model matrix will have ", N, " observations and ", P, " columns.\n")
        warning("too few observations to estimate. set poly_reg(..., warn_level=2) to stop on this error or poly_reg(..., warn_level=0) to bypass.")
      }
    }
    if(maxInteractDeg > 0)
      f <- paste(c(f, rep("*.", maxInteractDeg)), collapse="")

    return(model_matrix(formula(f), Xy))

  }

}
