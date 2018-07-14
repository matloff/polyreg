# polyreg

*Two Goals*

1.  Development of a package to automate formation and evaluation of
    polynomial regression models.  

2.  Comparison to neural networks, on the speculation that the latter
    operate like polynomial regression models.  Details in [Polynomial 
    Regression As an Alternative to Neural Nets](https://arxiv.org/abs/1806.06850), 
    by Cheng, Khomtchouk, Matloff and Mohanty, 2018

*Usage*

Other than the various cross-validation functions, the main functions
are **polyfit()** and **predict.polyFit()**.

Example: Programmer/engineer 2000 Census data, Silicon Valley.  Built in
to the latest version of [the **regtools**
package](https://github.com/matloff/regtools).  Predict wage income.

``` r
getPE()  # get dataset 
# try simple example, only a few predictors; wageinc lavt
pe <- pe[,c(1,2,4,6,7,3)]
pfout <- polyFit(pe,2)  # quadratic model
# predict wage of person age 40, male, 52 weeks worked, BS degree; need
# in data frame form, same names
newx <- pe[1,]  # dummy 1-row data frame
newx$age <- 40
newx$sex <- 1
newx$wkswrkd <- 52
newx <- newx[,-6]  # no Y value
predict(pfout,newx)  # about $68J
```
