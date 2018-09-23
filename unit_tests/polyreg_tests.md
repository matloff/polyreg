Unit Tests for polyreg
================
9/22/2018

This document tests out key features of `polyreg`.

Here are some tests for `FSR()`, the forward stepwise regression function. Under the hood, the `block_solve()` algorithm eases memory use and should yield the same results as `solve()` (up to rounding).

``` r
X <- as.matrix(mtcars)
baseR_approach  <- solve(crossprod(X)) 
baseR_approach[1:5, 1:5]
```

                   mpg           cyl          disp            hp          drat
    mpg   6.643240e-03 -0.0023306241 -8.996797e-05  1.365031e-04 -0.0082481281
    cyl  -2.330624e-03  0.0858854902 -7.072373e-04 -7.856340e-04  0.0007529713
    disp -8.996797e-05 -0.0007072373  4.660727e-05 -3.096288e-05 -0.0004012459
    hp    1.365031e-04 -0.0007856340 -3.096288e-05  6.998600e-05  0.0001312703
    drat -8.248128e-03  0.0007529713 -4.012459e-04  1.312703e-04  0.3229589439

``` r
polyreg_approach1 <- polyreg:::block_solve(X=X)
```

    ..

``` r
polyreg_approach1[1:5, 1:5]
```

                  [,1]          [,2]          [,3]          [,4]          [,5]
    [1,]  6.643240e-03 -0.0023306241 -8.996797e-05  1.365031e-04 -0.0082481281
    [2,] -2.330624e-03  0.0858854902 -7.072373e-04 -7.856340e-04  0.0007529713
    [3,] -8.996797e-05 -0.0007072373  4.660727e-05 -3.096288e-05 -0.0004012459
    [4,]  1.365031e-04 -0.0007856340 -3.096288e-05  6.998600e-05  0.0001312703
    [5,] -8.248128e-03  0.0007529713 -4.012459e-04  1.312703e-04  0.3229589439

``` r
max(abs(baseR_approach - polyreg_approach1)) < 10^{-11}
```

    [1] TRUE

``` r
polyreg_approach2 <- polyreg:::block_solve(S=crossprod(X))
```

    ..

``` r
polyreg_approach2[1:5, 1:5]
```

                  [,1]          [,2]          [,3]          [,4]          [,5]
    [1,]  6.643240e-03 -0.0023306241 -8.996797e-05  1.365031e-04 -0.0082481281
    [2,] -2.330624e-03  0.0858854902 -7.072373e-04 -7.856340e-04  0.0007529713
    [3,] -8.996797e-05 -0.0007072373  4.660727e-05 -3.096288e-05 -0.0004012459
    [4,]  1.365031e-04 -0.0007856340 -3.096288e-05  6.998600e-05  0.0001312703
    [5,] -8.248128e-03  0.0007529713 -4.012459e-04  1.312703e-04  0.3229589439

``` r
max(abs(baseR_approach - polyreg_approach2)) < 10^{-11}
```

    [1] TRUE

When `FSR()` estimates coefficients via Ordinary Least Squares, it should yield the same results as `lm()`.

``` r
library(polyreg)
baseR_beta <- coef(lm(carb ~ ., mtcars))
baseR_beta 
```

    (Intercept)         mpg         cyl        disp          hp        drat 
    -2.46807501 -0.01378803  0.28536857 -0.01431005  0.01349808  0.41696616 
             wt        qsec          vs          am        gear 
     1.53320915 -0.22493808 -0.23036244 -0.11878278  0.77153891 

``` r
X <- cbind(1, as.matrix(mtcars[,-ncol(mtcars)]))
y <- mtcars$carb
XtX_inv <- polyreg:::block_solve(X = X, max_block = 250)
```

    ..

``` r
polyreg_beta <- tcrossprod(XtX_inv, X) %*% y
polyreg_beta
```

                 [,1]
     [1,] -2.46807501
     [2,] -0.01378803
     [3,]  0.28536857
     [4,] -0.01431005
     [5,]  0.01349808
     [6,]  0.41696616
     [7,]  1.53320915
     [8,] -0.22493808
     [9,] -0.23036244
    [10,] -0.11878278
    [11,]  0.77153891

``` r
max(abs(baseR_beta - polyreg_beta)) < 10^{-9}
```

    [1] TRUE

`FSR()` works for different input types and uses different types of estimation. Specifically, continuous, binary, and multinomial outcomes may be estimated using OLS (as implemented above). Binary outcomes may also be estimated via `glm` and multinomial via `multinom::nnet`. Here is a simple way to make sure all five scenarios are running...

``` r
l <- FSR(mtcars, noisy = FALSE)
```

    .........

``` r
summary(l)
```

    The dependent variable is 'carb' which will be treated as continuous.  The data contains 32 observations (N_train == 24 and N_test == 8), which were split using seed 320174167. The data contains 8 continuous features and 0 dummy variables. Between 8 and 23 models will be estimated. Each model will add a feature, which will be included in subsequent models if it explains at least an additional 0.01 of variance out-of-sample (after adjusting for the additional term on [0, 1]).


    Estimated  8  models. 

    The best model has Predicted Adjusted R^2:
    The best model has Mean Absolute Predicted Error: 0.8546792

    The output has a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:

    predict(z, newdata = Xnew, model_to_use = 3) 

``` r
B <- cbind(mtcars, as.factor(sample(letters[1:2], nrow(mtcars), replace=TRUE)))
colnames(B)[ncol(B)] <- "idk"

b1 <- FSR(B, noisy=FALSE)
summary(b1)
```

    The dependent variable is 'idk' which will be treated as binary.  The data contains 32 observations (N_train == 27 and N_test == 5), which were split using seed 970843296. The data contains 9 continuous features and 1 dummy variables. Between 10 and 26 models will be estimated. Each model will add a feature, which will be included in subsequent models if it explains at least an additional 0.01 of variance out-of-sample (after adjusting for the additional term on [0, 1]).


    Estimated  10  models. 

    The best model has pseudo R^2 (adjusted for P and N):
    The best model has out-of-sample accuracy:

    The output has a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:

    predict(z, newdata = Xnew, model_to_use = 3) 

``` r
b2 <- FSR(B, linear_estimation = TRUE, noisy=FALSE)
```

    ...........

``` r
summary(b2)
```

    The dependent variable is 'idk' which will be treated as binary.  The data contains 32 observations (N_train == 24 and N_test == 8), which were split using seed 450151420. The data contains 9 continuous features and 1 dummy variables. Between 10 and 23 models will be estimated. Each model will add a feature, which will be included in subsequent models if it explains at least an additional 0.01 of variance out-of-sample (after adjusting for the additional term on [0, 1]).


    Estimated  10  models. 

    The best model has pseudo R^2 (adjusted for P and N):
    The best model has out-of-sample accuracy:

    The output has a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:

    predict(z, newdata = Xnew, model_to_use = 3) 

``` r
c1 <- FSR(iris, noisy=FALSE)
summary(c1)
```

    The dependent variable is 'Species' which will be treated as multinomial.  The data contains 150 observations (N_train == 121 and N_test == 29), which were split using seed 696000585. The data contains 4 continuous features and 2 dummy variables. Between 6 and 105 models will be estimated. Each model will add a feature, which will be included in subsequent models if it explains at least an additional 0.01 of variance out-of-sample (after adjusting for the additional term on [0, 1]).


    Estimated  6  models. 

    The best model has out-of-sample accuracy (adjusted for P and N): 0.9666667

    The output has a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:

    predict(z, newdata = Xnew, model_to_use = 3) 

``` r
c2 <- FSR(iris, linear_estimation = TRUE, noisy=FALSE)
```

    .......

``` r
summary(c2)
```

    The dependent variable is 'Species' which will be treated as multinomial.  The data contains 150 observations (N_train == 115 and N_test == 35), which were split using seed 514477689. The data contains 4 continuous features and 2 dummy variables. Between 6 and 105 models will be estimated. Each model will add a feature, which will be included in subsequent models if it explains at least an additional 0.01 of variance out-of-sample (after adjusting for the additional term on [0, 1]).


    Estimated  6  models. 

    The best model has out-of-sample accuracy (adjusted for P and N): 0.3964912

    The output has a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:

    predict(z, newdata = Xnew, model_to_use = 3)
