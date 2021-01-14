##################################################################
# getPoly: generate poly terms of a data matrix / data frame
##################################################################

#' @export
getPoly <- function(xdata = NULL, deg = 1, maxInteractDeg = deg,
                    Xy = NULL, 
                    standardize = FALSE,
                    noisy = TRUE, intercept = FALSE, 
                    returnDF = TRUE, 
                    modelFormula = NULL, retainedNames = NULL,
                    ...){

# maxInteractDeg <- maxInteractDeg + 1

  if(sum(is.null(xdata) + is.null(Xy)) != 1)
    stop("please provide getPoly() xdata or Xy (but not both).")
  

  W <- if(is.null(xdata)) Xy else xdata
  nX <- if(is.null(xdata)) ncol(Xy)-1 else ncol(xdata)  # NM, 11/13/2020
  if(noisy && !(is.matrix(W) || is.data.frame(W))){
    message("getPoly() expects a matrix or a data.frame. The input will be coerced to a data.frame but you may wish to stop and provide one directly.\n\n")
  }
  W <- as.data.frame(W, stringsAsFactors=TRUE)
  namesW <- names(W)
  W <- complete(W, noisy=noisy)
  colnames(W) <- namesW
  
  if(standardize){
    to_z <- which(unlist(lapply(W, is_continuous)))
    W[,to_z] <- scale(W[,to_z])
  }
  
  RHS <- NULL # initialize

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

      #if(W_distincts[i] == 2){
        tmp <- paste(colnames(W)[i], "==",
                     paste0("\'", levels(W[,i])[-1], "\'"))
        tmp <- paste0("(", tmp, ")")
        factor_features <- c(factor_features, tmp)
      #}else{
      #  factor_features <- c(factor_features, colnames(W)[i])
      #}
    }

    features <- factor_features
    P <- P_factor

    if(sum(x_factors) != ncol(W)){ # at least some continuous x variables
      
      # names of continuous variables
      cf <- colnames(W)[x_cols][!x_factors]
      P_continuous <- length(cf)
      P <- P_continuous + P_factor
      # P does not reflect intercept, interactions, or polynomial terms
      continuous_features <- c()
     for(i in 1:deg)
       continuous_features <- c(continuous_features, paste0("I(", cf, "^", i, ")"))
      # generate "I(x^1)", "I(x^2)", "I(x^3)", and so on
      # the string above will be used to make the appropriate formula
      # y ~ I(x^1) + I(x^2) + I(x^3)
      # "I(x^1)" is overkill but aids with debugging ... 
      features <- c(continuous_features, factor_features)
    } else cf <- NULL  # added by NM, 12/12/20
    
    if (ncol(W) > 1)
    features <- get_interactions(features, maxInteractDeg,
                                 c(cf, names(x_factors[x_factors])),
                                 maxDeg = deg)
    # get_interactions now returns original features too by default
    if(noisy && (length(features) > nrow(W)))
      message("P > N. With polynomial terms and interactions, P is ",
              length(features), ".\n\n", sep="")
    
    RHS <- paste0("~",  paste(features, collapse=" + "))
    modelFormula <- as.formula(paste0(y_name, RHS))
    # intercept handled by logical passed to model_matrix
    
  }
 
  ## if (ncol(W) == 1) names(W)[1] <- 'V1'  # bad kludge, NM, 11/12/20
  X <- model_matrix(modelFormula, W, intercept, noisy, ...)
  # if (!is.matrix(X)) X <- matrix(X,ncol=nX)
  if (is.vector(X)) {
     nms <- names(X)
     X <- matrix(X,nrow=nrow(W))
     colnames(X) <- nms
  }
  
  if(is.null(retainedNames))
    retainedNames <- colnames(X)
  
  polyMat <- polyMatrix(xdata = X, 
                        modelFormula = modelFormula, 
                        XtestFormula = if(is.null(RHS)) modelFormula else as.formula(RHS), 
                        retainedNames = retainedNames)

  if(returnDF) return(polyDF(polyMat)) else return(polyMat)

}

gtPoly <- getPoly

##################################################################
# polyMatrix: the class of polyMatrix from getPoly
##################################################################
# helper function so that getPoly() can also be called in parallel
# make polyMatrix <- function(...) so that it can be more flexible?
polyMatrix <- function(xdata, modelFormula, XtestFormula, retainedNames){
  
  if(!is.matrix(xdata)){
    xdata <- matrix(t(xdata), nrow=1)
    colnames(xdata) <- retainedNames
  }
  out <- list(xdata = xdata,
              modelFormula = modelFormula,
              XtestFormula = XtestFormula,
              retainedNames = retainedNames)
  class(out) <- "polyMatrix"
  return(out)

}

# helper function to coerce a polyMatrix to a polyDataFrame
polyDF <- function(polyMat){

  stopifnot(class(polyMat) == "polyMatrix")
  polyMat$xdata <- as.data.frame(polyMat$xdata, stringsAsFactors=polyMat$stringsAsFactors)
  if(length(polyMat$retainedNames) != length(colnames(polyMat$xdata))){
    warning("xdata contains", length(colnames(polyMat$xdata)), 
        "columns but retainedNames has", length(polyMat$retainedNames), "items")
  }
  class(polyMat) <- c("polyMatrix", "polyDataFrame")
  return(polyMat)
  
}


# parallel version of getPoly()
# needs updating
# or to be absorbed into model_matrix()
#
#getPolyPar <- function(cls, xdata)
#{
#  distribsplit(cls,'xdata')
#  cmd <- paste0('clusterEvalQ(cls, getPoly(xdata,',
#                deg,',',
##                maxInteractDeg,'))')
#  clusterEvalQ(cls,library(polyreg))
#  res <- eval(parse(text=cmd))
#  xd <- NULL
#  for (i in 1:length(cls)) {
#    xd <- rbind(xd,res[[i]]$xdata)
#  }
#  return (polyMatrix(xd, res[[1]]$endCols))
#}

polyAllVsAll <- function(plm.xy, classes, stringsAsFactors=TRUE){
  plm.xy <- as.data.frame(plm.xy, stringsAsFactors=TRUE)
  len <- length(classes)
  ft <- list()
  for (i in 1:len) {
    ft[[i]] <- list()
    for (j in 1:len) {
      if (i == j) # same class
        next
      newxy <- plm.xy[plm.xy$y == classes[i] | plm.xy$y == classes[j],]
      newxy$y <- ifelse(newxy$y == classes[i], 1, 0)
      ft[[i]][[j]] <- glm(y~., family = binomial(link = "logit"), data = newxy)
    } # for j
  } # for i
  return(ft)
}

polyOneVsAll <- function(plm.xy, classes, cls=NULL) {
  plm.xy <- as.data.frame(plm.xy, stringsAsFactors=TRUE)
  ft <- list()
  predClassi <- function(i)
  {
    oneclass <- plm.xy[plm.xy$y == classes[i],]
    oneclass$y <- 1
    allclass <- plm.xy[plm.xy$y != classes[i],]
    allclass$y <- 0
    new_xy <- rbind(oneclass, allclass)
    glm(y~., family = binomial(link = "logit"), data = new_xy)
  }
  if (is.null(cls)) {
    for (i in 1:length(classes)) {
      #     oneclass <- plm.xy[plm.xy$y == classes[i],]
      #     oneclass$y <- 1
      #     allclass <- plm.xy[plm.xy$y != classes[i],]
      #     allclass$y <- 0
      #     new_xy <- rbind(oneclass, allclass)
      #     ft[[i]] <- glm(y~., family = binomial(link = "logit"), data = new_xy)
      ft[[i]] <- predClassi(i)
    }
  } else {
    clusterExport(cls,c('plm.xy','predClassi'),envir=environment())
    ft <- clusterApply(cls,1:length(classes),predClassi)
  }
  return(ft)
}
