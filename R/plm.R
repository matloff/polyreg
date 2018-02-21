

##################################################################
# combnDeg: generate the degree distribution of different X's
##################################################################

# arguments:

#   n: number of X's
#   deg: total degrees

# return: a matrix of all possible degree distribution

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
    for (i in 2:col)
    {
      xy <- cbind(xy, as.numeric(unlist(xyd[i])))
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
# getPoly: generate poly terms of a data matrix / data frame
##################################################################

# arguments:

#   xydata: the dataframe (includes both predictor and response variables)
#   deg: the max degree of polynomial terms
#   maxInteractDeg: the max degree of dummy and nondummy predictor variable
#                   interaction terms

# return: a dataframe of all poly terms

#' @export
getPoly <- function(xydata, deg, maxInteractDeg = deg) {
  ### xydata includes y!

  if (deg < 1) {
    return("deg must be larger than or equal to 1.")

  }

  xydata <- as.data.frame(xydata)

  xy <- xydata[,-ncol(xydata), drop=FALSE]
  y <- xydata[,ncol(xydata)]
  n <- ncol(xy)
  # seperate dummy variables and continuous variables

  is_dummy <- (lapply(lapply(xy, table), length)==2)

  dummy <-xy[, is_dummy, drop = FALSE]
  nondummy <- xy[, !is_dummy, drop = FALSE]


  result <- xy


  if (deg > 1) {

    for (i in 2:deg) {
      if (ncol(nondummy) > 0) # for nondummy case
        result <- cbind(result, deg_plm(nondummy,i))


      if (ncol(dummy) > 0 && i <= ncol(dummy)){ # for dummy case
        result <- cbind(result, only_dummy(dummy,i))

      }
    }
  }

  if (maxInteractDeg > 1) {

    for (i in 2:maxInteractDeg) {
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
    } # for (i in 2:maxInteractDeg)
  }

  rt <- as.data.frame(cbind(result, y))
  colnames(rt)[ncol(rt)] <- "y"
  for (i in 1:(ncol(rt) - 1)) {
    colnames(rt)[i] <- paste("V", i, sep = "")
  }
  return (as.data.frame(rt))

}


##################################################################
# one_all: generate accuracy of classifier using one-vs-all method
##################################################################

# arguments:

#   xy: the data matrix/ dataframe
#   trainidxs: the indexes of training set

# return: the accuracy of the classification

one_all <- function(xy, trainidxs) {
  classes <- unique(xy[,ncol(xy)])
  training <- xy[trainidxs,]
  testing <- xy[-trainidxs,]
  prb <- matrix(0, nrow = nrow(testing),ncol=length(classes))
  for (i in 1:length(classes)) {
    oneClass <- training[training$y == classes[i],]
    oneClass$y <- 1
    allClass <- training[training$y != classes[i],]
    allClass$y <- 0
    new_xy <- as.data.frame(rbind(oneClass, allClass))
    pol <- glm(y~., family = binomial(link = "logit"), data = new_xy)
    pred <- predict(pol, testing[, -ncol(testing), drop=FALSE], na.action = na.omit, type = "response")
    prb[,i] <- pred
  }
  ridx <- apply(prb, 1, which.max)
  predF <- NULL
  for (k in 1:length(ridx)) {
    predF[k] <- classes[ridx[k]]
  }
  return(mean(predF == testing$y))
}


##################################################################
# all_all: generate accuracy of classifier using all-vs-all method
##################################################################

# arguments:

#   xy: the data matrix/ dataframe
#   trainidxs: the indexes of training set

# return: the accuracy of the classification

all_all <- function(xy, trainidxs) {
  classes <- unique(xy[,ncol(xy)])
  training <- xy[trainidxs,]
  testing <- xy[-trainidxs,]
  votes <- matrix(0, nrow = nrow(testing),ncol=length(classes))
  for (i in 1:length(classes)) {
    for (j in 1:length(classes)) {
      if (i == j) # same class
        next
      newtrain <- training[training$y == classes[i] | training$y == classes[j],]
      newtrain$y <- ifelse(newtrain$y == classes[i], 1, 0)
      pol <- glm(y~., family = binomial(link = "logit"), data = newtrain)
      pred <- predict(pol, testing[, -ncol(testing), drop=FALSE])
      predT <- ifelse(pred > 0.5, 1, 0)
      votes[,i] <- predT
    } # for j
  } # for i
  ridx <- apply(votes, 1, which.max)
  predF <- NULL
  for (k in 1:length(ridx)) {
    predF[k] <- classes[ridx[k]]
  }
  return(mean(predF == testing$y))
}


##################################################################
# polyFit: generate polynomial terms of data and fit models
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
polyFit <- function(xy, maxDeg, maxInteractDeg=maxDeg, use = "lm", trnProp=0.8,pcaMethod = FALSE,pcaPortion = 0.9, glmMethod = "all") {
  set.seed(500)
  n <- nrow(xy)
  ntrain <- round(trnProp*n)
  trainidxs <- sample(1:n, ntrain, replace = FALSE)

  error <- NULL
  for (i in 1:maxDeg) {  # for each degree
    m <- ifelse(i > maxInteractDeg, maxInteractDeg, i)
    pxy <- getPoly(xy, i, m)

    if (pcaMethod == TRUE) {
      pxy.pca <- prcomp(pxy[,-ncol(pxy)])
      pcNo = cumsum(pxy.pca$sdev)/sum(pxy.pca$sdev)
      for (k in 1:length(pcNo)) {
        if (pcNo[k] >= pcaPortion)
          break
      }
      trainingX <- pxy.pca$x[trainidxs,1:k]
      y <- pxy[trainidxs, ncol(pxy)]
      training <- as.data.frame(cbind(trainingX, y))
      testingX <- pxy.pca$x[-trainidxs,1:k]
      y <- pxy[-trainidxs, ncol(pxy)]
      testing <- as.data.frame(cbind(testingX, y))
      pxy <- as.data.frame(cbind(pxy.pca$x[,1:k], pxy[,ncol(pxy)]))
      colnames(pxy)[ncol(dt)] <- "y"
    } else if (pcaMethod == FALSE){
      training <- pxy[trainidxs,]
      testing <- pxy[-trainidxs,]
    } else {
      return("Please choose whether to use PCA (TRUE or FALSE).")
    }

    if (use == "lm") {
      pol <- lm(y~., data = training)
      pred <- predict(pol, testing[,-ncol(testing)], na.action = na.omit)
      error[i] <- mean(abs(pred - testing$y))
    }
    else if (use == "glm") {
      classes <- unique(xy[,ncol(xy)])
      if (length(classes) == 2) { # two categories
        training$y <- ifelse(training$y == classes[1], 1, 0)
        testing$y <- ifelse(testing$y == classes[1], 1, 0)
        pol <- glm(y~., family = binomial(link = "logit"), data = training)
        pred <- predict(pol, testing[, -ncol(testing)])
        predT <- ifelse(pred > 0.5, 1, 0)
        error[i] <- mean(predT == testing$y)
      }
      else { # more than two categories

        if (glmMethod == "one") {
          error[i] <- one_all(pxy, trainidxs)
        }
        else if (glmMethod == "all")
        {
          error[i] <- all_all(pxy, trainidxs)
        }
        else {
          return("Please choose 'one' for one-vs-all or 'all' for all-vs-all. ")
        }
      } # else more than two categories

    } # glm
    else {
      return("Please choose lm or glm.")

    }
  }
  return(error)

}




