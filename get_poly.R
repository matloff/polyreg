##################################################################
# get_poly: generate poly terms of a data matrix / data frame
##################################################################

# arguments:

#   xdata: the dataframe (only predictor variables). Factors with more than two levels should not be inputted as integers.
#   deg: the max degree of polynomial terms.
#   maxInteractDeg: the max degree of dummy and nondummy predictor variable
#      interaction terms
#   Xy: the dataframe with the response in the final column (provide xdata or Xy but not both).
#     Factors with more than two levels should not be inputted as integers.
#   model_formula: Internal use. Formula used to generate the training model matrix.
#     Note: anticipates that polynomial terms are generated using internal functions of library(polyreg)
#     so YMMV if the formula is not generated on the training data by get_poly().
#     Also, providing model_formula bypasses deg and maxInteractDeg.
#   standardize: standardize all continuous variables? (Default: FALSE.)
#   intercept: Include intercept? Default: FALSE.
# return: a model matrix, with the model formula as an additional attribute
# examples:
# X = get_poly(mtcars, 2)
# W = get_poly(Xy=mtcars, 2) # treats final column as response
# ncol(W) < ncol(X)          # TRUE
# X_train <- get_poly(mtcars[1:20,], 4, 2)
# X_test <- get_poly(mtcars[21:32,],
#                    model_formula = attributes(X_train)$formula)
get_poly <- function(xdata = NULL, deg=1, maxInteractDeg = deg,
                     Xy = NULL, model_formula = NULL, standardize = FALSE,
                     intercept = FALSE){

  if(sum(is.null(xdata) + is.null(Xy)) != 1)
    stop("please provide get_poly() xdata or Xy (but not both).")

  W <- if(is.null(xdata)) as.data.frame(Xy) else as.data.frame(xdata)
  if(standardize){
    to_z <- which(unlist(lapply(W, is_continuous)))
    W[,to_z] <- scale(W[,to_z])
  }


  if(is.null(model_formula)){

    x_cols <- if(is.null(xdata)) 1:(ncol(Xy) - 1) else 1:ncol(xdata)
    y_name <- if(is.null(xdata)) colnames(Xy)[ncol(Xy)] else NULL

    remove(xdata)
    remove(Xy)

    # coerce binary or character variables into factors
    W_distincts <- N_distinct(W)
    to_factor <- which((W_distincts == 2) | unlist(lapply(W, is.character)))
    for(i in to_factor)
      W[,i] <- as.factor(W[,i])
    x_factors <- unlist(lapply(W[,x_cols], is.factor))
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

    continuous_features <- cf <- colnames(W)[x_cols][!x_factors]
    P_continuous <- length(continuous_features)

    P <- P_continuous + P_factor
    # P does not reflect intercept, interactions, or polynomial terms

    for(i in 2:deg)
      continuous_features <- c(continuous_features,
                               paste("pow(", cf, ",", i, ")"))
    # pow() is a helper function that deals with the nuissance
    # that lm(y ~ x + x^2 + x^3)
    # will only estimate one slope but we want three...
    # the string above will be used to make the appropriate formula
    # y ~ x + pow(x, 2) + pow(x, 3)

    features <- c(continuous_features, factor_features)

    if(maxInteractDeg > 1)
      features <- get_interactions(features, maxInteractDeg, names(x_factors[x_factors]))

    if(length(features) > nrow(W))
      warning("P > N. With polynomial terms and interactions, P is ",
              length(features))

    model_formula <- as.formula(paste0(y_name, " ~ ",
                                       ifelse(intercept, "", "-1 +"),
                                      paste(features, collapse=" + ")))

  }

  return(model_matrix(model_formula, W))

}
