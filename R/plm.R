
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
# endCols:  a list endCols[i] would be the column number of xy in
#             which the degree-i terms end
#' @export
polyMatrix <- function(x, k) {
  me <- list(xdata = x, endCols = k)
  class(me) <- "polyMatrix"
  return(me)
}



##################################################################
# getPoly: generate poly terms of a data matrix / data frame
##################################################################

# arguments:

#   xdata: the dataframe (only predictor variables)
#   deg: the max degree of polynomial terms
#   maxInteractDeg: the max degree of dummy and nondummy predictor variable
#                   interaction terms

# return: a polyMatrix object

#' @export
getPoly <- function(xdata, deg, maxInteractDeg = deg)
{

  if (deg < 1) {
    return("deg must be larger than or equal to 1.")
  }

  xdata <- as.data.frame(xdata)
  endCols <- NULL
  #xy <- xydata[,-ncol(xydata), drop=FALSE]
  #y <- xydata[,ncol(xydata)]
  n <- ncol(xdata)
  # seperate dummy variables and continuous variables
  is_dummy <- (lapply(lapply(xdata, table), length)==2)
  dummy <-xdata[, is_dummy, drop = FALSE]
  nondummy <- xdata[, !is_dummy, drop = FALSE]

  result <- xdata
  endCols[1] <- ncol(xdata) # deg 1 starts at result[1] (first column)

  if (deg > 1) {

    for (m in 2:deg) {
      i <- ifelse(m > maxInteractDeg, maxInteractDeg, m)

      if (ncol(nondummy) > 0) # for nondummy case
        result <- cbind(result, deg_plm(nondummy,m))

      if (ncol(dummy) > 0 && i <= ncol(dummy)) # for dummy case
        result <- cbind(result, only_dummy(dummy,m))

      # for dummy & nondummy intersection
      if (ncol(nondummy) > 0 && ncol(dummy) > 0) {
        for (j in 1:(i-1)) {
          if (j == 1 && i - j == 1) {
            r_dummy <- dummy
            r_nondummy <- nondummy
          }
          else if (j == 1) { # when dummy is only distributed 1 deg
            r_dummy <- dummy
            r_nondummy <- deg_plm(nondummy,i-j)
          }
          else if (i - j == 1) { # the case when nondummy is
                                 # only distributed 1 deg
            r_dummy <- only_dummy(dummy, j)
            r_nondummy <- nondummy
          }
          else {
            r_nondummy <- deg_plm(nondummy,i-j)
            r_dummy <- only_dummy(dummy, j)
          }
          mix <- as.numeric(r_dummy[,1]) * as.numeric(r_nondummy[,1])
          skip <- 1
          n_dummy <- ncol(r_dummy)
          n_nondummy <- ncol(r_nondummy)
          for (a in 1:n_dummy) {
            for (b in 1:n_nondummy) {
              if (skip == 1) {
                skip <- skip - 1
                next
              }
              mix <- cbind(mix,
                 as.numeric(r_dummy[,a]) * as.numeric(r_nondummy[,b]))
            }
          }
          result <- cbind(result, mix)
        }
      } # dummy & nondummy intersection

      endCols[m] <- ncol(result)
    } # loop 2:deg
  } # if deg > 1

  rt <- as.data.frame(result)
  for (i in 1:ncol(rt)) {
    colnames(rt)[i] <- paste("V", i, sep = "")
  }

  return (polyMatrix(rt, endCols))
}


polyAllvsAll <- function(plm.xy, classes){
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

polyOnevsAll <- function(plm.xy, classes) {
  plm.xy <- as.data.frame(plm.xy)
  ft <- list()
  for (i in 1:length(classes)) {
    oneClass <- plm.xy[plm.xy$y == classes[i],]
    oneClass$y <- 1
    allClass <- plm.xy[plm.xy$y != classes[i],]
    allClass$y <- 0
    new_xy <- rbind(oneClass, allClass)
    ft[[i]] <- glm(y~., family = binomial(link = "logit"), data = new_xy)
  }
  return(ft)
}

##################################################################
# polyFit: generate polynomial terms of data and fit models
##################################################################

# arguments:
#   xy: dataframe, response variable is in the last column
#   deg: the degree of the polynomial terms
#   maxInteractDeg: the max degree of dummy and nondummy predictor variables
#                   interaction terms
#   use: can be "lm" for linear regreesion, and "glm" for logistic regression
#   pcaMethod: whether to use PCA (can be TRUE or FALSE)
#   pcaPortion: the portion of principal components to be used
#   glmMethod: which method ("all" for all-vs-all, "one" for one-vs-all)
#              to use for multi-class classification
#   printTimes: whether to print the time of PCA, getPoly, lm, or glm.
#   polyMat: if non-NULL, then polynomial matrix will be passed in

# return: the object of class polyFit

#' @export
polyFit <- function(xy, deg, maxInteractDeg, use = "lm", pcaMethod = FALSE,
     pcaPortion = 0.9, glmMethod = "all",printTimes=TRUE, polyMat = NULL) {
  y <- xy[,ncol(xy)]
  classes <- FALSE
  xy.pca <- NULL
  k <- 0
  if (!is.null(polyMat)) { # polynomial matrix is provided
    if (pcaMethod == TRUE) {
      tmp <- system.time(
        xy.pca <- prcomp(polyMat[,-ncol(xy)])
      )
      if (printTimes) cat('PCA time: ',tmp,'\n')
      pcNo = cumsum(xy.pca$sdev)/sum(xy.pca$sdev)
      for (k in 1:length(pcNo)) {
        if (pcNo[k] >= pcaPortion)
          break
      }
      xdata <- xy.pca$x[,1:k, drop=FALSE]
      
    } else {
      xdata <- polyMat
    }
    plm.xy <- cbind(xdata,y)
  } # polynomial matrix is provided
  else { # polynomial matrix is not provided
    if (pcaMethod == TRUE) {
      tmp <- system.time(
        xy.pca <- prcomp(xy[,-ncol(xy)])
      )
      if (printTimes) cat('PCA time: ',tmp,'\n')
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
      plm.xy <- cbind(getPoly(xdata, deg, maxInteractDeg)$xdata,y)
    )
    if (printTimes) cat('getPoly time: ',tmp,'\n')
  } # polynomial matrix is not provided
  

  if (use == "lm") {
    tmp <- system.time(
    ft <- lm(y~., data = plm.xy)
    )
    if (printTimes) cat('lm() time: ',tmp,'\n')
    glmMethod <- NULL
  }
  else if (use == "glm") {
    classes <- unique(y)
    if (length(classes) == 2) {
      plm.xy$y <- as.numeric(ifelse(plm.xy$y == classes[1], 1, 0))
      tmp <- system.time(
      ft <- glm(y~., family = binomial(link = "logit"), data = plm.xy)
      )
      if (printTimes) cat('2-class glm() time: ',tmp,'\n')
      glmMethod <- NULL
    }
    else { # more than two classes
      if (glmMethod == "all") { # all-vs-all
        tmp <- system.time(
        ft <- polyAllvsAll(plm.xy, classes)
        )
        if (printTimes) cat('all-vs-all glm() time: ',tmp,'\n')
      } else if (glmMethod == "one") { # one-vs-all
        tmp <- system.time(
        ft <- polyOnevsAll(plm.xy, classes)
        )
        if (printTimes) cat('one-vs-all glm() time: ',tmp,'\n')
      }
    } # more than two classes

  }
  pcaPrn <- ifelse(pcaMethod == TRUE, pcaPortion, 0)
  me <- list(xy = xy, degree = deg, maxInteractDeg = maxInteractDeg, use = use,
             poly.xy = plm.xy, fit = ft, PCA = pcaMethod, pca.portion = pcaPrn,
             pca.xy = xy.pca, pcaCol = k, glmMethod = glmMethod,
             classes = classes)
  class(me) <- "polyFit"
  return(me)

}

##################################################################
# predict.polyFit: predict the fitted models on newdata
##################################################################

# arguments:
#   object: a polyFit class object
#   newdata: a new dataframe, without response variable
#   polyMat: if non-NULL, then polynomial matrix will be passed in

# return: predicted values of newdata

#' @export
predict.polyFit <- function(object, newdata, polyMat = NULL) { 
  # newdata doesn't have y column
  
  if (!is.null(polyMat)) { # polynomial matrix is provided
    if (object$PCA == TRUE) {
      plm.newdata <- predict(object$pca.xy, polyMat)[,1:object$pcaCol]
    } else {
      plm.newdata <- polyMat
    }
    plm.newdata <- as.data.frame(plm.newdata)
  } 
  else { # polynomial matrix is not provided
    if (object$PCA == TRUE) {
      new_data <- predict(object$pca.xy, newdata)[,1:object$pcaCol]
      plm.newdata <- getPoly(new_data, object$degree, object$maxInteractDeg)$xdata
    } else {
      plm.newdata <- getPoly(newdata, object$degree, object$maxInteractDeg)$xdata
    }
  } # polynomial matrix is not provided
  

  if (object$use == "lm") {
    pred <- predict(object$fit, plm.newdata)
  } else { # glm case
    if (is.null(object$glmMethod)) { # only two classes
      pre <- predict(object$fit, plm.newdata)
      pred <- ifelse(pre > 0.5, object$classes[1], object$classes[2])
    } else { # more than two classes
      len <- length(object$classes)
      if (object$glmMethod == "all") { # all-vs-all method
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
          prob[,i] <- predict(object$fit[[i]], plm.newdata, type = "response")
        }
        winner <- apply(prob, 1, which.max)
      }
      pred <- NULL
      for (k in 1:nrow(plm.newdata)) {
        pred[k] <- object$classes[winner[k]]
      }
    } # more than two classes
  } # glm case

  return(pred)
}

