

##################################################################
# combnDeg: generate the degree distribution of different X's
##################################################################

# arguments:

#   n: number of X's
#   deg: total degree

# return: a matrix of all possible degree distributions

# say only have X1, X2, i.e. 2 predictors and 'deg' degree;
# then the degree distribution can first be: 1, deg-1, i.e. the 
# exponents of X1^i * X2^j, and similarly
#  2, deg-2
#  3, deg-3
#  ...
#  deg-1, 1
combnDeg <- function(n, deg) { # distribute (deg) degrees to (n) different X's

  if (n == 2) { 
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
      mydata[[i]] <- matrix(i, nrow=nrow(temp), ncol=1) # set the degree for the first variable
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
# startCols:  a list startCols[i] would be the column number of xy in
#             which the degree-i terms begin
#' @export
polyMatrix <- function(x, k) {
  me <- list(xy = x, endCols = k)
  class(me) <- "polyMatrix"
  return(me)
}



##################################################################
# getPoly: generate poly terms of a data matrix / data frame
##################################################################

# arguments:

#   xdf: the dataframe (only predictor variables)
#   deg: the max degree of polynomial terms
#   maxInteractDeg: the max degree of dummy and nondummy predictor variable
#                   interaction terms

# return: a polyMatrix object

#' @export
getPoly <- function(xdf, deg, maxInteractDeg = deg) {
  ### xdf doesn't include y

  if (deg < 1) {
    return("deg must be larger than or equal to 1.")
  }

  xdf <- as.data.frame(xdf)
  endCols <- NULL
  #xy <- xydata[,-ncol(xydata), drop=FALSE]
  #y <- xydata[,ncol(xydata)]
  n <- ncol(xdf)
  # seperate dummy variables and continuous variables
  is_dummy <- (lapply(lapply(xdf, table), length)==2)
  dummy <-xdf[, is_dummy, drop = FALSE]
  nondummy <- xdf[, !is_dummy, drop = FALSE]

  result <- xdf
  endCols[1] <- ncol(xdf) # deg 1 starts at result[1] (first column)

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
          else if (i - j == 1) { # the case when nondummy is only distributed 1 deg
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
              mix <- cbind(mix, as.numeric(r_dummy[,a]) * as.numeric(r_nondummy[,b]))
            }
          }
          result <- cbind(result, mix)
        }
      } # dummy & nondummy intersection

      endCols[m] <- ncol(result)
    } # loop 2:deg
  } # if deg > 1

  rt <- as.data.frame(result)
  #colnames(rt)[ncol(rt)] <- "y"
  for (i in 1:ncol(rt)) {
    colnames(rt)[i] <- paste("V", i, sep = "")
  }

  return (polyMatrix(rt, endCols))

}


polyAllvsAll <- function(plm.xy, classes){
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

# return: the object of class polyFit, S3 object with components:

# xy,degree,maxInteractDeg,use,PCA(pcaMethod),pca.portion(pcaPortion,glmMethod:
#    the input arguments of the same/similar names
# poly.xy:  cbind(poly matrix of the Xs + Y vector)
# fit:  output of lm() or glm(), except in multiclass case
# pca.xy:  output of prcomp() on the Xs matrix
# pcaCol:  number of principal components selected according to pcaPortion

#' @export
polyFit <- function(xy, deg, maxInteractDeg, use = "lm", pcaMethod=FALSE, 
     pcaPortion = 0.9, glmMethod = "all") 
{
  y <- xy[,ncol(xy)]
  classes <- FALSE
  rotation <- FALSE
  xy.pca <- NULL
  k <- 0
  if (pcaMethod == TRUE) {
    xy.pca <- prcomp(xy[,-ncol(xy)])
    pcNo = cumsum(xy.pca$sdev)/sum(xy.pca$sdev)
    for (k in 1:length(pcNo)) {
      if (pcNo[k] >= pcaPortion)
        break
    }
    xdata <- xy.pca$x[,1:k, drop=FALSE]

  } else {
    xdata <- xy[,-ncol(xy), drop=FALSE]
  }
  plm.xy <- cbind(getPoly(xdata, deg, maxInteractDeg)$xy,y)

  if (use == "lm") {
    ft <- lm(y~., data = plm.xy)
    glmMethod <- NULL
  }
  else if (use == "glm") {
    classes <- unique(y)
    if (length(classes) == 2) {
      plm.xy$y <- as.numeric(ifelse(plm.xy$y == classes[1], 1, 0))
      ft <- glm(y~., family = binomial(link = "logit"), data = plm.xy)
      glmMethod <- NULL
    }
    else { # more than two classes
      if (glmMethod == "all") { # all-vs-all
        ft <- polyAllvsAll(plm.xy, classes)
      } else if (glmMethod == "one") { # one-vs-all
        ft <- polyOnevsAll(plm.xy, classes)
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

# return: predicted values of newdata

#' @export
predict.polyFit <- function(object, newdata) { # newdata doesn't have y column
  if (object$PCA == TRUE) {
    # newdata <- (scale(newdata,scale=FALSE, center=TRUE) 
    # %*% object$pcaRotation)[,1:object$pcaCol]
    new_data <- predict(object$pca.xy, newdata)[,1:object$pcaCol]
    plm.newdata <- getPoly(new_data, object$degree, object$maxInteractDeg)$xy
  } else {
    plm.newdata <- getPoly(newdata, object$degree, object$maxInteractDeg)$xy
  }

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
xvalPoly <- function(xy, maxDeg, maxInteractDeg=maxDeg, 
      use = "lm", trnProp=0.8,pcaMethod = FALSE,pcaPortion = 0.9, 
      glmMethod = "all") 
{
  set.seed(500)
  n <- nrow(xy)
  ntrain <- round(trnProp*n)
  trainidxs <- sample(1:n, ntrain, replace = FALSE)
  error <- NULL
  training <- xy[trainidxs,]
  testing <- xy[-trainidxs,]

  for (i in 1:maxDeg) {  # for each degree
    m <- ifelse(i > maxInteractDeg, maxInteractDeg, i)
    pol <- polyFit(training, i, m, use, pcaMethod, pcaPortion, glmMethod)
    pred <- predict(pol, testing[,-ncol(testing)])
    if (use == "lm") {
      error[i] <- mean(abs(pred - testing[,ncol(testing)])) # abs. mean error
    } else if (use == "glm")
      error[i] <- mean(pred == testing[,ncol(testing)]) # accuracy

  }
  return(error)

}

