##################################################################
# polyFit: generate polynomial terms of data and fit models
##################################################################

# arguments:
#   xy: dataframe, response variable is in the last column; in
#      classification case (indicated by 'use', this must be either 
#      an R factor or a numeric code for the various classes
#   deg: the degree of the polynomial terms
#   maxInteractDeg: the max degree of dummy and nondummy predictor variables
#      interaction terms
#   use: can be "lm" for linear regreesion, "glm" for logistic
#      regression, or "mvrlm" for multivariate-response lm(); the
#      latter two indicate a classification case
#   glmMethod: which method ("all" for all-vs-all, "one" for one-vs-all,
#      "multlog" for multinomial logistic regression)
#      to use for multi-class classification
#   cls:  R 'parallel' cluster; currently not used

# return: the object of class polyFit

# note on RSpectra(): if the full set of eigenvectors is requested,
# RSpectra() simply calls R's built-in eigen(), so we don't treat that
# case here

# NM, 09/19/18: removed dropout option, getting in the way and not
# useful

polyFit <- function(xy, deg, maxInteractDeg=deg, use = "lm", 
                    glmMethod="one", return_xy = FALSE, 
                    returnPoly = FALSE, noisy = TRUE)
{

  if (!use %in% c('lm','glm','mvrlm'))
     stop('"use" must be "lm", "glm", or "mvrlm"')
 
  nOrigFeatures <- ncol(xy) - 1
  namesOrigFeatures <- colnames(xy[,1:nOrigFeatures,drop=FALSE])

  xy <- complete(xy, noisy=noisy)
  
  xdata <- xy[,-ncol(xy),drop=FALSE]

  y <- xy[,ncol(xy)]
  if(is.character(y))
    y <- as.factor(y)
  # is this a classification problem?
  classProblem <- is.factor(y) || use == 'mvrlm'
  if (classProblem) {
     classes <- levels(y)
     if (is.factor(y))  { # change to numeric code for the classes
        y <- as.numeric(y)
        numClasses <- 1:length(classes)
        xy[,ncol(xy)] <- y
     }
  } else classes <- FALSE

     k <- 0
     tmp <- system.time(pMat <- getPoly(xdata, deg, maxInteractDeg))
     if(noisy) message('getPoly time: ', max(tmp, na.rm = TRUE),'\n\n')
     polyMat <- pMat$xdata

  retainedNames <- pMat$retainedNames
  modelFormula <- pMat$modelFormula
  XtestFormula <- pMat$XtestFormula

  # by now, polyMat is ready for input to lm() etc. in all cases

  # this is the new xy, i.e. the polynomialized and possibly PCA-ized
  # version of xy
  plm.xy <- as.data.frame(cbind(polyMat,y), stringsAsFactors=TRUE)

  # OK, PCA and getPoly() taken care of, now find the fit, to be
  # assigned to ft

  if (use == "lm") {
    tmp <- system.time(
       ft <- lm(y~., data = plm.xy)
    )
    if(noisy) message('lm() time: ', max(tmp, na.rm=TRUE),'\n\n')
    glmMethod <- NULL
  } else if (use == "glm" || use == 'mvrlm') {
       if (use == 'glm') {
          if (length(numClasses) == 2) {
            plm.xy$y <- as.numeric(plm.xy$y == numClasses[1])
            tmp <- system.time(ft <- glm(y~., family = binomial,data = plm.xy))
            if(noisy) 
              message('2-class glm() time: ', max(tmp, na.rm = TRUE),'\n\n')
            glmMethod <- NULL
          }  # end 2-class case
          else { # more than two classes
            if (glmMethod == "all") { # all-vs-all
              tmp <- system.time(ft <- polyAllVsAll(plm.xy, numClasses))
              if(noisy) 
                message('all-vs-all glm() time: ', max(tmp,na.rm=TRUE),'\n\n')
            } else if (glmMethod == "one") { # one-vs-all
              tmp <- system.time(
                 ft <- polyOneVsAll(plm.xy, numClasses) # cls could be passed here
              )
              glmOuts <- ft
              if(noisy) 
                 message('one-vs-all glm() time: ', max(tmp,na.rm=TRUE),'\n\n')
            } else if (glmMethod == "multlog") { # multinomial logistics
              tmp <- system.time(
              ft <- multinom(y~., plm.xy)
              )
              if(noisy) message('multlog time: ', max(tmp, na.rm = TRUE), '\n')
            }
          } # more than two classes
      # end 'glm' case
      }  else  {  # 'mvrlm' case
            # dms <- model.matrix(~ as.factor(y) - 1, y)
            yf <- as.factor(y)
            dms <- model.matrix(~yf-1)
            dms <- as.data.frame(dms, stringsAsFactors=TRUE)
            dxy <- cbind(plm.xy[,-ncol(plm.xy)],dms)
            nms <- names(dms)
            addnames <- paste0(nms,collapse=',')
            frml <- paste0('cbind(',addnames,') ~ .,data=dxy')
            # somehow as.formula() has a problem here, so back to basics
            cmd <- paste0('ft <- lm(',frml,')')
            eval(parse(text=cmd))
      }

  }  # end 'glm'/'mvrlm' case

  if (!exists('glmOuts')) glmOuts <- NULL
  
  me <- list(xy = if(return_xy) xy else NULL, 
             degree=deg, 
             maxInteractDeg=maxInteractDeg, 
             use=use,
             poly.xy = if(returnPoly) plm.xy else NULL, 
             fit=ft, 
             nOrigFeatures=nOrigFeatures,
             namesOrigFeatures=namesOrigFeatures,
             glmMethod=glmMethod,
             classProblem=classProblem, 
             classes=classes, 
             retainedNames=retainedNames, 
             modelFormula=modelFormula, 
             XtestFormula=XtestFormula,
             glmOuts=glmOuts)
  class(me) <- "polyFit"
  return(me)

}



