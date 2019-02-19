


#### 09/11/18, NM:  overhauling this file, kind of a mess

##################################################################
# combnDeg: generate the degree distribution of different X's
##################################################################

# arguments:

#   n: number of X's
#   deg: total degrees

# return: a matrix of all possible degree distribution

# example1: if n = 2, deg = 3 -> 2 predictor variables and 3 degrees
#          this function would return:
#                 [ 1, 2,
#                   2, 1 ]
# for interaction terms: X1, X2, and maxDeg=3
# then use the first row as the first possible degrees: X1*X2^2
# also use the second row as the second possible degrees: X1^2*X2
# Therefore, 2 predictors & 3 degrees would have interaction terms:
#            X1 * X2^2 + X1^2 * X2

# example2: if n = 3, deg = 4 -> 3 predictor variables and 4 degrees
#           this function would return:
#                  [ 1, 1, 2,
#                    1, 2, 1,
#                    2, 1, 1 ]
# coresponding to 3 interaction terms:
#  X1 * X2 * X3^2 + X1 * X2^2 * X3 + X1^2 * X2 * X3

# further examples in the comments at the head of getPoly()

combnDeg <- function(n, deg) { # distribute (deg) degrees to (n) different X's

  if (n == 2) { # if only have X1, X2, *2* predictors and *deg* degrees
    # then the degree distribution can be: 1, deg-1
    #  2, deg-2
    #  3, deg-3
    #  ...
    #  deg-1, 1
    result <- matrix(0, nrow=deg-1, ncol=2)
    for (i in 1:(deg-1)) {
      result[i,] <- c(i, deg-i)
    }
    return (result)
  } # if n==2 (base case)
  else if (n > 2) { # if have more than 2 predictors
    # set the degree for the first variable, recursive on the
    # rest eg. X1,X2,X3, deg=4 then X1 can have deg=1, and
    # X2,X3 can have a total of deg=3 (only two variables, go
    # to base case)      X1 can have deg=2, and X2,X3 can have
    # a total of deg=2      ...
    mydata <- list()
    for (i in 1:(deg-n+1)) {
      temp <- combnDeg(n-1,deg-i)
      # set the degree for the first variable
      mydata[[i]] <- matrix(i, nrow=nrow(temp), ncol=1)
      mydata[[i]] <- cbind(mydata[[i]], temp)
    }
    result <- mydata[[1]]
    # do the first one seperately because it needs the first one to cbind later

    if (length(mydata) > 1) {
      for (j in 2:(deg-n+1)) {
        result <- rbind(result, mydata[[j]])
        # combine the rows of different degrees for first variable
      }
    }
    return (result)
  }
  else {
    print("Error on combnDeg.")
  }
}


##################################################################
# deg_plm: generate poly terms for nondummy variables
##################################################################

# arguments:

#   xyd: the data matrix of nondummy predictor variables
#   deg: the degree of polynomial terms

# return: the new data matrix with polynomial terms

deg_plm <- function(xyd, deg) { # deal with nondummy terms only

  row <- dim(xyd)[1]
  col <- dim(xyd)[2]

  if (typeof(xyd) == "list") {
    xy <- data.frame(as.numeric(unlist(xyd[1])))
    if (col > 1) {
      for (i in 2:col)
      {
        xy <- cbind(xy, as.numeric(unlist(xyd[i])))
      }
    }
  } # fix type
  else {
    xy <- xyd
  }
  result <- xy^deg


  lim <- min(col, deg)
  if (col > 1) {
    for (i in 2:lim) {
      idx <- combn(1:col, i) # get the combination of col index
      # (i.e. different combination of predictors)
      idx_row <- nrow(idx)
      idx_col <- ncol(idx)

      if (i <= deg)
        deg_dist <- combnDeg(i, deg)
      else
        deg_dist <- combnDeg(deg, deg)
      # get the different distributions of degrees on each predictor
      deg_row <- nrow(deg_dist)
      deg_col <- ncol(deg_dist)

      for (j in 1:deg_row) {
        for(k in 1:idx_col) {
          temp <- matrix(1, ncol = 1, nrow = row)
          for (l in 1:idx_row) {
            # choose variable using col index
            # and choose its degree using the distribution of degrees
            temp <- temp * xy[, idx[l,k] ]^(deg_dist[j,l])
          }
          result <- cbind(result, temp)
        }
      }
    }
  }

  return (result)
}



##################################################################
# only_dummy: generate poly terms for dummy variables
##################################################################

# arguments:

#   xy: the data matrix of dummy predictor variables
#   deg: the degree of polynomial terms

# return: the new data matrix with polynomial terms

only_dummy <- function(xy, deg) { # deal with dummy terms only

  n <- ncol(xy) # number of predictors

  if (n <= deg) { # deal with the case: eg. X1,X2, deg=3
    # --> X1*X2^2 = X1^2*X2 = just need X1*X2
    result <- matrix(1, ncol=1, nrow=nrow(xy))
    for (i in 1:n) {
      result <- result * as.numeric(xy[,i])
    }
  }
  else { # if n > deg, deal with the case: eg. X1,X2,X3, deg=2
    # --> choose two of them
    idx <- combn(1:n, deg) # get different combinations of variables

    result <- matrix(1,ncol=1, nrow=nrow(xy))
    for (k in 1:nrow(idx)) {
      result <- result * as.numeric(xy[, idx[k,1]])
    } # do the first one seperately
    # because it needs the first one to cbind later

    if (ncol(idx) ==1)
      return(result)

    for (j in 2:ncol(idx)) {
      temp <- matrix(1,ncol=1, nrow=nrow(xy))
      for (k in 1:nrow(idx))
        temp <- temp * as.numeric(xy[, idx[k,j]])

      result <- cbind(result, temp)
    }
  }

  return(result)
}

##################################################################
# polyMatrix: the class of polyMatrix from getPoly
##################################################################

# xy:  the actual xy matrix
# zeroDel:  all-0 columns created during polynomial generation were deleted
polyMatrix <- function(x, retainedNames) {
  me <- list(xdata = x, retainedNames = retainedNames)
  class(me) <- "polyMatrix"
  return(me)
}


# deletes all-0 columns of the data frame d, returning the new data frame
delete0columns <- function(d)
{
   all0 <- function(x) all(x == 0)
   tmp <- apply(d,2,all0)
   colsKeep <- which(!tmp)
   d[,colsKeep]
}

# parallel version of getPoly()
getPolyPar <- function(cls,xdata,deg,maxInteractDeg=deg)
{
   distribsplit(cls,'xdata')
   cmd <- paste0('clusterEvalQ(cls,getPoly(xdata,',
             deg,',',
             maxInteractDeg,'))')
   # clusterExport(cls,c('deg_plm','combnDeg','getPoly'),envir=environment())
   clusterEvalQ(cls,library(polyreg))
   res <- eval(parse(text=cmd))
   xd <- NULL
   for (i in 1:length(cls)) {
      xd <- rbind(xd,res[[i]]$xdata)
   }
   return (polyMatrix(xd, res[[1]]$endCols))
}

polyAllVsAll <- function(plm.xy, classes){
  plm.xy <- as.data.frame(plm.xy)
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

polyOneVsAll <- function(plm.xy, classes,cls=NULL) {
  plm.xy <- as.data.frame(plm.xy)
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

# to test whenever changes are made to getPoly(); output should be

#    a c bu bv
# 1  4 2  1  0
# 2 14 3  0  1
# 3 17 4  0  0
# 4 10 5  1  0
# 5 12 6  0  0
# > gpx <- getPoly(x,2)
# > gpx
# $xdata
#   V1 V2 V3 V4  V5 V6 V7 V8 V9 V10 V11 V12
# 1  4  2  1  0  16  4  8  0  4   2   0   0
# 2 14  3  0  1 196  9 42  0  0   0  14   3
# 3 17  4  0  0 289 16 68  0  0   0   0   0
# 4 10  5  1  0 100 25 50  0 10   5   0   0
# 5 12  6  0  0 144 36 72  0  0   0   0   0
#
# $endCols
# [1]  4 12
#
# attr(,"class")
# [1] "polyMatrix"

testGP <- function()
{
#   requireNamespace(dummies)
   x <-
      data.frame(a=sample(1:20,5),b=factor(c('u','v','w','u','w')),c=2:6)
   tmp <- dummy(x$b)
   x <- cbind(x,tmp[,-3])
   x$b <- NULL
}

##################################################################
# polyFit: generate polynomial terms of data and fit models
##################################################################

# arguments:
#   xy: dataframe, response variable is in the last column; in
#      classification case, this must be either an R factor or a
#      numeric code for the various classes
#   deg: the degree of the polynomial terms
#   maxInteractDeg: the max degree of dummy and nondummy predictor variables
#      interaction terms
#   use: can be "lm" for linear regreesion, "glm" for logistic
#      regression, or "mvrlm" for multivariate-response lm()
#   pcaMethod: default is NULL, can be either "prcomp" (use the prcomp()
#      function to compute PCA) or "RSectra" (use sparse Matrix and
#      compute eigenvalues/vectors to compute PCA)
#   pcaPortion: number of principal components to be used; if < 1, this
#      specifies a desired proportion of explained variance,
#      otherwise the actual number of components; if RSectra
#      method is used, must be >= 1
#   pcaLocation: if 'front', compute principal comps and then form
#      polynomial in them; if 'back', do the opposite;
#      relevant only if pcaMethod is non-NULL
#   glmMethod: which method ("all" for all-vs-all, "one" for one-vs-all,
#      "multlog" for multinomial logistic regression)
#      to use for multi-class classification
#   cls:  R 'parallel' cluster

# return: the object of class polyFit

# note on RSpectra(): if the full set of eigenvectors is requested,
# RSpectra() simply calls R's built-in eigen(), so we don't treat that
# case here

# NM, 09/19/18: removed dropout option, getting in the way and not
# useful

polyFit <- function(xy, deg, maxInteractDeg=deg, use = "lm", pcaMethod=NULL,
     pcaLocation='front', pcaPortion=0.9, glmMethod="one", cls=NULL)
{

  if (!use %in% c('lm','glm','mvrlm'))
     stop('"use" must be "lm", "glm" or "mvrlm"')

  doPCA <- !is.null(pcaMethod)
  xdata <- xy[,-ncol(xy)]

  y <- xy[,ncol(xy)]
  if(is.character(y))
    y <- as.factor(y)
  # is this a classification problem?
  classProblem <- is.factor(y) || use == 'mvrlm'
  if (classProblem) {
     if (is.factor(y))  { # change to numeric code for the classes
        y <- as.numeric(y)
        xy[,ncol(xy)] <- y
     }
     classes <- unique(y)
  } else classes <- FALSE

  if (doPCA)  {  # start PCA section
    # safety checks first
    if (pcaMethod == 'RSpectra' && pcaPortion < 1)
       stop('use prcomp method for this case')
    if (!pcaMethod %in% c('prcomp','RSpectra'))
       stop("pcaMethod should be either NULL, prcomp, or RSpectra")
    stopifnot(pcaLocation %in% c('front','back'))
    # can't do PCA with R factors or char
    if (!all(apply(xdata,2,is.numeric)))
       stop('X data must be numeric for PCA')
    # now compute
    if (pcaLocation == 'front') {
       applyPCAOutputs <- applyPCA(xdata,pcaMethod,pcaPortion)
       xdata <- applyPCAOutputs$xdata
       tmp <-
         system.time(pMat <- getPoly(xdata, deg, maxInteractDeg))
       cat('getPoly time: ',tmp,'\n')
  polyMat <- pMat$xdata
  retainedNames <- pMat$retainedNames
    } else  {  # 'back'

      tmp <- system.time(pMat <- getPoly(xdata, deg, maxInteractDeg))
      cat('getPoly time: ',tmp,'\n')
      polyMat <- pMat$xdata
      retainedNames <- pMat$retainedNames
      applyPCAOutputs <- applyPCA(polyMat,pcaMethod,pcaPortion)
      polyMat <- applyPCAOutputs$xdata

    }
    xy.pca <- applyPCAOutputs$xy.pca  # overall output of prcomp or RSpectra
    k <- applyPCAOutputs$k
  # end PCA section
  }  else {   # no-PCA section
     xy.pca <- NULL
     k <- 0
     tmp <- system.time(pMat <- getPoly(xdata, deg, maxInteractDeg))
     cat('getPoly time: ',tmp,'\n')
  polyMat <- pMat$xdata
  retainedNames <- pMat$retainedNames
  }


  # by now, polyMat is ready for input to lm() etc. in all cases

  # this is the new xy, i.e. the polynomialized and possibly PCA-ized
  # version of xy
  plm.xy <- as.data.frame(cbind(polyMat,y))

  # OK, PCA and getPoly() taken care of, now find the fit, to be
  # assigned to ft

  if (use == "lm") {
    tmp <- system.time(
       ft <- lm(y~., data = plm.xy)
    )
    cat('lm() time: ',tmp,'\n')
    glmMethod <- NULL
  } else if (use == "glm" || use == 'mvrlm') {
       classes <- unique(y)  # see preprocessing of y, start of this ftn
       if (use == 'glm') {
          if (length(classes) == 2) {
            # plm.xy$y <- as.numeric(ifelse(plm.xy$y == classes[1], 1, 0))
            plm.xy$y <- as.numeric(plm.xy$y == classes[1])
            tmp <- system.time(ft <- glm(y~., family = binomial,data = plm.xy))
            cat('2-class glm() time: ',tmp,'\n')
            glmMethod <- NULL
          }  # end 2-class case
          else { # more than two classes
            if (glmMethod == "all") { # all-vs-all
              tmp <- system.time(ft <- polyAllVsAll(plm.xy, classes))
              cat('all-vs-all glm() time: ',tmp,'\n')
            } else if (glmMethod == "one") { # one-vs-all
              tmp <- system.time(
                 ft <- polyOneVsAll(plm.xy, classes,cls)
              )
              cat('one-vs-all glm() time: ',tmp,'\n')
            } else if (glmMethod == "multlog") { # multinomial logistics
               #requireNamespace(nnet)
              tmp <- system.time(
              ft <- multinom(y~., plm.xy)
              )
              cat('multlog time: ', tmp, '\n')
            }
          } # more than two classes
      # end 'glm' case
      }  else  {  # 'mvrlm' case
            #requireNamespace(dummies)
            dms <- dummy(y)
            dms <- as.data.frame(dms)
            dxy <- cbind(plm.xy[,-ncol(plm.xy)],dms)
            nms <- names(dms)
            addnames <- paste0(nms,collapse=',')
            frml <- paste0('cbind(',addnames,') ~ .,data=dxy')
            # somehow as.formula() has a problem here, so back to basics
            cmd <- paste0('ft <- lm(',frml,')')
            eval(parse(text=cmd))
      }

  }  # end 'glm'/'mvrlm' case

  # create return value and wrap up
  pcaPrn <- if(doPCA) pcaPortion else 0
  me <-list(xy=xy,degree=deg,maxInteractDeg=maxInteractDeg,use=use,
    poly.xy=plm.xy,fit=ft,PCA=pcaMethod,pca.portion=pcaPrn,
    pca.xy=xy.pca,pcaCol=k,pcaLocation=pcaLocation,glmMethod=glmMethod,
    classProblem=classProblem,classes=classes,retainedNames=retainedNames)
  class(me) <- "polyFit"
  return(me)

}

# 09/11/18, NM: moved this function out of polyFit(), now standalone,
# for readability

applyPCA <- function(x,pcaMethod,pcaPortion) {

  if (pcaMethod == "prcomp") { # use prcomp for pca
    tmp <- system.time(
      #xy.pca <- prcomp(x[,-ncol(xy)])
      xy.pca <- prcomp(x)
    )
    cat('PCA time: ',tmp,'\n')
    if (pcaPortion >= 1.0) k <- pcaPortion else {
       k <- 0
       pcNo = cumsum(xy.pca$sdev)/sum(xy.pca$sdev)
       for (k in 1:length(pcNo)) {
         if (pcNo[k] >= pcaPortion)
           break
       }
    }
    cat(k,' principal comps used\n')
    xdata <- xy.pca$x[,1:k, drop=FALSE]

  } else { # use RSpectra for PCA
    #requireNamespace(RSpectra)
    xy.cov <- cov(x)
    k <- pcaPortion
    xy.eig <- eigs(xy.cov,k)
    xy.pca <- xy.eig
    cat(k,' principal comps used\n')
    #xdata <- as.matrix(x[,-ncol(x)]) %*% xy.eig$vectors[,1:k]
    xdata <- as.matrix(x) %*% xy.eig$vectors[,1:k]
  }

  return(list(xdata=xdata,xy.pca=xy.pca,k=k))
}

##################################################################
# predict.polyFit: predict the fitted models on newdata
##################################################################

# NM, 09/17/18: removed the polyMat argument; unlikely to be used,
# cluttered up the code

# arguments:
#   object: a polyFit class object
#   newdata: a new dataframe, same column names as in the training
#            data but without response variable

# return: predicted values of newdata, IN THE FORM OF NUMERICAL CLASS
# CODES, 1,2,3,...
predict.polyFit <- function(object, newdata, ...)
{
  use <- object$use

  # the next couple dozen lines are devoted to forming plm.newdata, which
  # will ultimately be fed into predict.lm(), predict.glm() or whatever;
  # to do this, newdata, the argument above, must be expanded to
  # polynomial form, and possibly run through PCA

  doPCA <- !is.null(object$PCA)

  if (!doPCA) {

    browser()

    f <- paste("~", paste(object$retainedNames, collapse=" + "))
    plm.newdata <- getPoly(newdata, object$degree, object$maxInteractDeg, modelFormula = f, ...)$xdata

  } else if (object$PCA == "prcomp") {

        message("Beginning PCA\n\n", timestamp())

       if (object$pcaLocation == "front") {

         new_data <- predict(object$pca.xy, newdata)[,1:object$pcaCol]
         plm.newdata <- getPoly(new_data, object$degree, object$maxInteractDeg, ...)$xdata

       } else if (object$pcaLocation == "back") {

         new_data <- getPoly(newdata, object$degree, object$maxInteractDeg, ...)$xdata
         plm.newdata <- predict(object$pca.xy, new_data)[,1:object$pcaCol]
         plm.newdata <- as.data.frame(plm.newdata)

       } else stop('invalid pcaLocation. Should be "front" or "back".')
  } else if (object$PCA == "RSpectra") {

       if (object$pcaLocation == "front") {

         xy.eig <- object$pca.xy
         new_data <- as.matrix(newdata) %*% xy.eig$vectors
         plm.newdata <- getPoly(new_data, object$degree, object$maxInteractDeg)$xdata

       } else if (object$pcaLocation == "back") {
         new_data <- getPoly(newdata, object$degree, object$maxInteractDeg, ...)$xdata
         ### xy.cov <- cov(new_data)
         ### xy.eig <- eigs(xy.cov,object$pcaCol)
         xy.eig <- object$pca.xy
         plm.newdata <-
            as.matrix(new_data) %*% xy.eig$vectors[,1:object$pcaCol]
         plm.newdata <- as.data.frame(plm.newdata)
       }
  }  # end doPCA
  message("Finished with PCA and model matrix construction.\n\n", timestamp())


  if (object$use == "lm") {
    pred <- predict(object$fit, plm.newdata)
    return(pred)
  }

  if (object$use == "mvrlm") {
      pre <- predict(object$fit, plm.newdata)
      pred <- apply(pre,1,which.max)
      return(pred)
  }

  # glm case
  if (is.null(object$glmMethod)) { # only two classes
    pre <- predict(object$fit, plm.newdata)
    pred <- ifelse(pre > 0.5, object$classes[1], object$classes[2])
  } else { # more than two classes
    len <- length(object$classes)
    if (object$glmMethod == "multlog") { # multinomial logistics
      pr <- predict(object$fit, plm.newdata, type="probs")
      idx <- apply(pr,1, which.max)
      col.name <- colnames(pr)
      lc <- length(col.name)
      tempM <- matrix(rep(col.name, length(idx)), ncol=lc, byrow = TRUE)
      pred <- NULL
      for (r in 1:nrow(tempM)) {
        pred[r] <- tempM[r,idx[r]]
      }
      return(pred)
    } # end multinomial logistics
    else if (object$glmMethod == "all") { # all-vs-all method
      votes <- matrix(0, nrow = nrow(plm.newdata), ncol = len)
      for (i in 1:len) {
        for (j in 1:len) {
          if (i == j)
            next
          pre <- predict(object$fit[[i]][[j]], plm.newdata, type="response")
          votes[,i] <- votes[,i] + ifelse(pre > 0.5, 1, 0)
        } # for i
      } # for j
      winner <- apply(votes, 1, which.max)

    } else if (object$glmMethod == "one") { # one-vs-all method
        prob <- matrix(0, nrow=nrow(plm.newdata), ncol=len)
        for (i in 1:len) {
          # prob[,i] <- parSapplyLB(object$fit[[i]],
          prob[,i] <- predict(object$fit[[i]],
             plm.newdata, type = "response")
        }
      winner <- apply(prob, 1, which.max)
    } # one-vs-all method
    # calculate pred for all-vs-all & one-vs-all
    pred <- NULL
    for (k in 1:nrow(plm.newdata)) {
      pred[k] <- object$classes[winner[k]]
    }
  } # end more than two classes
  return(pred)
  # end glm case

}

