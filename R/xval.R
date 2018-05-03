
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
                     yCol = NULL,printTimes=TRUE) 
{
  if (!is.null(yCol)) xy <- moveY(xy,yCol)
  
  y <- xy[,ncol(xy)]
  
  if (pcaMethod) {
    tmp <- system.time(
      xy.pca <- prcomp(xy[,-ncol(xy)])
    )
    if (printTimes) cat('PCA time in xvalPoly: ',tmp,'\n')
    pcNo = cumsum(xy.pca$sdev)/sum(xy.pca$sdev)
    for (k in 1:length(pcNo)) {
      if (pcNo[k] >= pcaPortion)
        break
    }
    xdata <- xy.pca$x[,1:k, drop=FALSE]
  } else {
    xdata <- xy[,-ncol(xy), drop=FALSE]
  }
  
  tmp <- system.time(
    poly.xy <- getPoly(xdata, maxDeg, maxInteractDeg)
  )
  if (printTimes) cat('getPoly time in xvalPoly: ',tmp,'\n')
  
  xy <- cbind(poly.xy$xdata, y)
  
  tmp <- splitData(xy,nHoldout)
  training <- tmp$trainSet
  testing <- tmp$testSet
  train.y <- training[,ncol(training)]
  test.y <- testing[,ncol(testing)]
#  testIdx <- splitData(xy, nHoldout,TRUE)
  
  acc <- NULL
  for (i in 1:maxDeg) {  # for each degree
    m <- ifelse(i > maxInteractDeg, maxInteractDeg, i)
    
    endCol <- poly.xy$endCols[i]
    
#    if (use == "glm" & glmMethod == "multlog" & length(unique(y)) > 2) {
#      # multinomial logistic
#      dat <- cbind(xy[,1:endCol], y)
#      pol <- polyFit(dat, i, m, use, pcaMethod = FALSE, pcaPortion, glmMethod,
#                     polyMat = dat, testIdx = testIdx)
#      pred <- predict(pol, test1, test1, testIdx)
#    } else {
      train1 <- cbind(training[,1:endCol], train.y)
      colnames(train1)[ncol(train1)] <- "y"
      test1 <- testing[,1:endCol]
      
      pol <- polyFit(train1, i, m, use, pcaMethod = FALSE, pcaPortion, glmMethod,
                     polyMat = train1)
      pred <- predict(pol, test1, test1)
#    }
    
    if (use == "lm") {
      acc[i] <- mean(abs(pred - test.y))
    } else
      acc[i] <- mean(pred == test.y) # accuracy
    
  } # for each degree
  return(acc)
}

##################################################################
# xvalNnet: generate mean absolute error of fitted models
##################################################################

# xval for nnet package

# arguments and value: 
#    see xvalPoly() above for most
#    classification is done if the Y variable is a factor
#    scaleXMat: if TRUE, apply scale() to predictor matrix
#    xy: as above, except that in classification case,
#        Y column of xy (last one, or yCol) must be a factor

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
  cmd <- paste0('nnout <- nnet(',yName,' ~ .,data=training,size=')
  cmd <- paste0(cmd,size,',linout=',linout,')')
  eval(parse(text=cmd))
  preds <- predict(nnout,testingx)
  trainingy <- training[,ncxy]
  if (!is.factor(trainingy))  # regression case
    return(mean(abs(preds - testingy)))
  # classification case
  preds <- apply(preds,1,which.max)  # column numbers
  # convert to levels of Y
  preds <- levels(trainingy)[preds]
  return(mean(preds == testingy))
}

##################################################################
# xvalKf: generate mean absolute error of fitted models
##################################################################

# xval for nnet package

# arguments and value: 
#    xy: as above, except that in classification case,
#        Y column of xy (last one, or yCol) must be a factor; otherwise
#        must NOT be a factor
#    nHoldout,yCol as above
#    scale_coninuous: "zero-one" for the classification case, otherwise NULL
#    loss:  last-stage activation ftn, '"mean_squared_error"' for
#           regression case, otherwise NULL
#    rmArgs:  remaining optional arguments to kms()

# return: a vector of mean absolute error (for lm) or accuracy (for glm),
#         the i-th element of the list is for degree = i
#' @export

xvalKf <- function(xy,nHoldout=10000,yCol=NULL,rmArgs=NULL)
{
  require(kerasformula)
  ncxy <- ncol(xy)
  
  if (!is.null(yCol)) xy <- moveY(xy,yCol)
  
  tmp <- splitData(xy,nHoldout)
  training <- tmp$trainSet
  testing <- tmp$testSet
  testingx <- tmp$testSet[,-ncxy]
  testingy <- tmp$testSet[,ncxy]
  
  
  yName <- names(xy)[ncol(xy)]
  trainingy <- training[,ncxy]
  classcase <- is.factor(trainingy)
  loss <- if (classcase) 'NULL' else "mean_squared_error"
  cmd <- paste0('kfout <- kms(',yName,' ~ .,data=training,loss=')
  cmd <- paste0(cmd,loss,')')
  eval(parse(text=cmd))
  preds <- predict(kfout,testingx)$fit
  if (!classcase) {  # regression case
     ry <- range(trainingy)
     preds <- ry[1] + (ry[2]-ry[1]) * preds
     return(mean(abs(preds - testingy)))
  } 
  # classification case
  return(mean(preds == testingy))
}

##################################################################
# xvalDNet: generate mean absolute error of fitted models
##################################################################

# xval for deepnet package, nn.*()

# arguments: 
#    not all args of nn.train() are implemented here, e.g. dropout args
#       are missing
#    hidden: vector of number of units in each layer
#    output: final layer feeds into this; double quoted;
#            '"sig"m', '"linear"' or '"softmax"'
#    numepochs: number of epochs
#    pca*, scaleXMat,nHoldout,yCol: as above
#    x: predictor variables
#    y: Y; in classification case, must be a matrix of dummies

# value: mean abs. error

#' @export

xvalDnet <- function(x,y,hidden,output='"sigm"',numepochs=3,
                     pcaMethod = FALSE,pcaPortion = 0.9,
                     scaleXMat = TRUE, nHoldout=10000)
{
  require(deepnet)
  
  if (scaleXMat) x <- scale(x)
  
  tmp <- splitData(x,nHoldout,idxsOnly=TRUE)
  trainingx <- x[-tmp,]
  testingx <- x[tmp,]
  ym <- as.matrix(y)
  trainingy <- ym[-tmp,]
  testingy <- ym[tmp,]
  
  cmd <- paste0('nnout <- nn.train(trainingx,trainingy,')
  cmd <- paste0(cmd,'hidden=',hidden,',') 
  cmd <- paste0(cmd,'output=',output,',') 
  cmd <- paste0(cmd,'numepochs=',numepochs) 
  cmd <- paste0(cmd,')')
  eval(parse(text=cmd))
  preds <- nn.predict(nnout,testingx)
  if (ncol(ym) == 1)  # regression case
    return(mean(abs(preds - testingy)))
  # else
  preds <- apply(preds,1,which.max)  # column numbers
  trueY <- apply(testingy,1,function(rw) which(rw == 1))
  # convert to levels of Y
  return(mean(preds == trueY))
  
}

######################  splitData() #################################
# support function, to split into training and test sets
##################################################################

splitData <- function(xy,nHoldout,idxsOnly=FALSE) 
{
  n <- nrow(xy)
  set.seed(500)
  ntrain <- nHoldout
  testIdxs <- sample(1:n, ntrain, replace = FALSE)
  if (idxsOnly) return(testIdxs)
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

# xy consists of X columns followed by one Y column, maybe a factor,
# unless xOnly is TRUE, in which case xy is just the X columns
scaleX <- function(xy) 
{  
  ncxy <- ncol(xy)
  x <- xy[,-ncxy]
  x <- scale(x)
  xy[,-ncxy] <- x
  xy
}
