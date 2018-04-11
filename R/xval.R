
# TODO:

#    xval*() should take advantage of the endCols component in the
#    polyMatrix class, in extending a lower-degree model to a
#    higher-degree one

##################################################################
# xvalPoly: generate mean absolute error of fitted models
##################################################################

# arguments:
#   xy: dataframe, response variable is in the last column
#   maxDeg: the max degree of the polynomial terms
#   maxInteractDeg: the max degree of dummy and nondummy predictor variables
#                   interaction terms
#   use: can be "lm" for linear regreesion, and "glm" for logistic regression
#   trnProp: the portion of data to be used as training set
#   pcaMethod: whether to use PCA (can be TRUE or FALSE)
#   pcaPortion: the portion of principal components to be used
#   glmMethod: when the response variable has more than 2 classes, can be "all"
#              for All-vs-All method, and "one" for One-vs-All method

# return: a list of mean absolute error (for lm) or accuracy (for glm),
#         the i-th element of the list is for degree = i
#' @export
xvalPoly <- function(xy, maxDeg, maxInteractDeg = maxDeg, use = "lm",
                     pcaMethod = FALSE,pcaPortion = 0.9, 
                     glmMethod = "all",nHoldout=10000,
                     yCol = NULL) 
{
  if (!is.null(yCol)) xy <- moveY(xy,yCol)
  tmp <- splitData(xy,nHoldout)
  training <- tmp$trainSet
  testing <- tmp$testSet

  error <- NULL
  for (i in 1:maxDeg) {  # for each degree
    m <- ifelse(i > maxInteractDeg, maxInteractDeg, i)
    pol <- polyFit(training, i, m, use, pcaMethod, pcaPortion, glmMethod)
    pred <- predict(pol, testing[,-ncol(testing)])
    if (use == "lm") {
      # absolute mean error
      error[i] <- mean(abs(pred - testing[,ncol(testing)]))
    } else if (use == "glm")
      error[i] <- mean(pred == testing[,ncol(testing)]) # accuracy

  }
  return(error)
}

##################################################################
# xvalNNet: generate mean absolute error of fitted models
##################################################################

######################  splitData() #################################
# support function, to split into training and test sets
##################################################################

splitData <- function(xy,nHoldout) 
{
  n <- nrow(xy)
  set.seed(500)
  ntrain <- nHoldout
  testIdxs <- sample(1:n, ntrain, replace = FALSE)
  testSet <- xy[testIdxs,]
  trainSet <- xy[-testIdxs,]
  list(testSet=testSet,trainSet=trainSet)
}

######################  moveY() #################################
# support function, since getPoly() etc. require Y in last column
#################################################################

moveY <- function(xy,yCol) 
{
  yName <- names(xy)[yCol]
  xy <- cbind(xy[,-yCol],xy[,yCol])
  names(xy)[ncol(xy)] <- yName
  xy
}

######################  scaleX() #################################
# support function, since NNs tend to like scaling
##################################################################

# scaleX <- function() 
# {
# }
