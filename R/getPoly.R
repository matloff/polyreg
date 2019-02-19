##################################################################
# getPoly: generate poly terms of a data matrix / data frame
##################################################################

#' @export
getPoly <- function(xdata = NULL, deg = 1, maxInteractDeg = deg,
                     Xy = NULL, modelFormula = NULL, standardize = FALSE,
                     noisy = TRUE, intercept = FALSE, ...){

  if(sum(is.null(xdata) + is.null(Xy)) != 1)
    stop("please provide getPoly() xdata or Xy (but not both).")

  W <- if(is.null(xdata)) as.data.frame(Xy) else as.data.frame(xdata)
  if(standardize){
    to_z <- which(unlist(lapply(W, is_continuous)))
    W[,to_z] <- scale(W[,to_z])
  }

  if(is.null(modelFormula)){

    x_cols <- 1:(ncol(W) - is.null(xdata))
    y_name <- if(is.null(xdata)) colnames(Xy)[ncol(Xy)] else NULL

    remove(xdata)
    remove(Xy)

    # coerce binary or character variables into factors
    W_distincts <- N_distinct(W)
    to_factor <- which((W_distincts == 2) | unlist(lapply(W, is.character)))
    #x_factors <- vector("logical", length(x_cols))
    for(i in to_factor)
      W[,i] <- as.factor(W[,i])
    x_factors <- if(ncol(W) > 1) unlist(lapply(W[,x_cols], is.factor)) else is.factor(W[[1]])
    P_factor <- sum(x_factors)

    factor_features <- c()
    # stores individual levels of factor variables, omitting one as reference
    # e.g. for a variable "sex" coded as binary, just "male" will be stored
    # this enables the appropriate formula to be written that handles interactions
    # suppose the variable party had three levels the strings
    # 'party == GOP' and 'party == independent' are stored
    # and democrat is the reference...

    for(i in which(x_factors)){

      if(W_distincts[i] > 2){
        tmp <- paste(colnames(W)[i], "==",
                     paste0("\'", levels(W[,i])[-1], "\'"))
        tmp <- paste0("(", tmp, ")")
        factor_features <- c(factor_features, tmp)
      }else{
        factor_features <- c(factor_features, colnames(W)[i])
      }
    }

    continuous_features <- cf <- colnames(W)[x_cols][!x_factors] # cf needed as copy below
    P_continuous <- length(continuous_features)

    P <- P_continuous + P_factor
    # P does not reflect intercept, interactions, or polynomial terms

    if(length(continuous_features) > 0){
      for(i in 2:deg)
        continuous_features <- c(continuous_features, paste0("I(", cf, "^", i, ")"))
    }
    # generate "I(x^2)", "I(x^3)", and so on
    # the string above will be used to make the appropriate formula
    # y ~ x + I(x^2) + I(x^3)

    features <- c(continuous_features, factor_features)

    #if(maxInteractDeg > 1 && ncol(W) > 1) # checks now in get_interactions
    features <- get_interactions(features, maxInteractDeg,
                                 c(cf, names(x_factors[x_factors])),
                                 maxDeg = deg)
    # get_interactions now returns original features too by default

    if(noisy && (length(features) > nrow(W)))
      cat("P > N. With polynomial terms and interactions, P is ",
              length(features), ".\n\n", sep="")

    modelFormula <- as.formula(paste0(y_name, " ~ ",
                                      # ifelse(intercept, "", "-1 +"), # https://stats.stackexchange.com/questions/174976/why-does-the-intercept-column-in-model-matrix-replace-the-first-factor
                                      paste(features, collapse=" + ")))

  }
  X <- model_matrix(modelFormula, W, intercept, noisy, ...)
  #if(!("formula" %in% names(attributes)))
  #  attr(X, which="formula") <- modelFormula # patch, should be addressed in model_matrix()
  return(list(xdata = X, retainedNames = colnames(X)))

}
