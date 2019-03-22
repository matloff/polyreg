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
       tmp <- system.time(pMat <- getPoly(xdata, deg, maxInteractDeg))
       cat('getPoly time: ',tmp,'\n')
       polyMat <- pMat$xdata
    #   retainedNames <- pMat$retainedNames
    } else  {  # 'back'

      tmp <- system.time(pMat <- getPoly(xdata, deg, maxInteractDeg))
      cat('getPoly time: ',tmp,'\n')
      polyMat <- pMat$xdata
   #   retainedNames <- pMat$retainedNames
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
  #   retainedNames <- pMat$retainedNames
  }
  retainedNames <- pMat$retainedNames
  modelFormula <- pMat$modelFormula
  XtestFormula <- pMat$XtestFormula

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
            # requireNamespace(dummies)
            # require(dummies)
            # dms <- dummies::dummy(y)
            # dms <- model.matrix(~ as.factor(y) - 1, y)
            yf <- as.factor(y)
            dms <- model.matrix(~yf-1)
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
  
  me <- list(xy=xy, degree=deg, maxInteractDeg=maxInteractDeg, use=use,
    poly.xy=plm.xy, fit=ft, PCA=pcaMethod, pca.portion=pcaPrn,
    pca.xy=xy.pca, pcaCol=k, pcaLocation=pcaLocation, glmMethod=glmMethod,
    classProblem=classProblem, classes=classes, 
    retainedNames=retainedNames, # retainedNames should perhaps be depracated with the new getPoly()
    modelFormula=modelFormula, 
    XtestFormula=XtestFormula)
  class(me) <- "polyFit"
  return(me)

}



