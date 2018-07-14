# polyreg, an Alternative to Neural Networks

Development of a package to automate formation and evaluation of
multivariate polynomial regression models.  

Motivation:  A simpler, equally effective alternative to neural
networks.  in [Polynomial Regression As an Alternative to Neural 
Nets](https://arxiv.org/abs/1806.06850), by Cheng, Khomtchouk, 
Matloff and Mohanty, 2018

**Usage:**

Other than the various cross-validation functions, the main functions
are **polyfit()** and **predict.polyFit()**.  One can fit either
regression or classification models, with an option to perform PCA for
dimension reduction on the predictors/features.

*Example:* Programmer/engineer 2000 Census data, Silicon Valley.  

Built in to the latest version of [the **regtools**
package](https://github.com/matloff/regtools).  Install package or
download directly 
[here](https://github.com/matloff/regtools/raw/master/data/prgeng.txt).
In the former case, **getPE()** reads in the dataset and does some
preprocessing, producing a data frame **pe**.

``` r
getPE()  # get dataset 
# predict wage income
# try simple example, only a few predictors; wageinc last
pe <- pe[,c(1,2,4,6,7,3)]
# take a look
head(pe,2)
#        age sex wkswrkd ms phd wageinc
# 1 50.30082   0      52  0   0   75000
# 2 41.10139   1      20  0   0   12300
pfout <- polyFit(pe,2)  # quadratic model
# predict wage of person age 40, male, 52 weeks worked, BS degree
# need in data frame form, same names
newx <- pe[1,]  # dummy 1-row data frame
newx <- newx[,-6]  # no Y value
newx$age <- 40
newx$sex <- 1
newx
#   age sex wkswrkd ms phd
# 1  40   1      52  0   0
predict(pfout,newx)  # about $68K
```
