
##################################################################
# xvalPoly: generate mean absolute error of fitted models
##################################################################

# arguments:
#   xy: dataframe, response variable is in the last column
#   maxDeg: the max degree of the polynomial terms
#   maxInteractDeg: the max degree of dummy and nondummy predictor variables
#                   interaction terms
#   use: can be "lm" for linear regreesion, and "glm" for logistic regression
#   pcaMethod: whether to use PCA (can be TRUE or FALSE)
#   pcaPortion: the portion of principal components to be used
#   glmMethod: when the response variable has more than 2 classes, can be "all"
#              for All-vs-All method, and "one" for One-vs-All method
#   nHoldout: number of cases for the test set
#   yCol: if not NULL, Y is in this column, and will be moved to last

# return: a vector of mean absolute error (for lm) or accuracy (for glm),
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

# arguments and value: 
#    see xvalPoly() above for most
#    scaleXMat: if TRUE, apply scale() to predictor matrix

# return: a vector of mean absolute error (for lm) or accuracy (for glm),
#         the i-th element of the list is for degree = i
#' @export

xvalNnet <- function(xy,size,linout, pcaMethod = FALSE,pcaPortion = 0.9,
               scaleXMat = FALSE, nHoldout=10000, yCol = NULL)
{
   require(nnet)
   ncxy <- ncol(xy)

   if (!is.null(yCol)) xy <- moveY(xy,yCol)
 
   if (scaleXMat) xy <- scaleX(xy)  # only Xs are scaled
   tmp <- splitData(xy,nHoldout)
   training <- tmp$trainSet
   testingx <- tmp$testSet[,-ncxy]
   testingy <- tmp$testSet[,ncxy]

   yName <- names(xy)[ncol(xy)]
   cmd <- paste0('nnout <- nnet(',yName,' ~ .,data=xy,size=')
   cmd <- paste0(cmd,size,',linout=',linout,')')
   eval(parse(text=cmd))
   npred <- predict(nnout,testingx)
   mean(abs(npred - testingy))
}

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

scaleX <- function(xy) 
{
   ncxy <- ncol(xy)
   nms <- names(xy)
   x <- xy[,-ncxy]
   x <- scale(x)
   tmp <- as.data.frame(cbind(x,xy[,ncxy]))
   names(tmp) <- nms
   tmp
}

