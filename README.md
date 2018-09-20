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

*Example:* Vertebral Column data from the [UC Irvine Machine
Learning Repository](https://archive.ics.uci.edu/ml/datasets/Vertebral+Column).
Various spinal measurements, with three conditions, Normal, 
Disk Hernia and Spondylolisthesis.  Let's predict the conditions.

``` r
# vert <- read.table('~/Research/DataSets/Vertebrae/column_3C.dat',header=FALSE)
vert$V7 <- as.character(vert$V7)  # Y must be a vector, not a factor
head(vert)
#      V1    V2    V3    V4     V5    V6 V7
# 1 63.03 22.55 39.61 40.48  98.67 -0.25 DH
# 2 39.06 10.06 25.02 29.00 114.41  4.56 DH
# 3 68.83 22.22 50.09 46.61 105.99 -3.53 DH
# 4 69.30 24.65 44.31 44.64 101.87 11.21 DH
# 5 49.71  9.65 28.32 40.06 108.17  7.92 DH
# 6 40.25 13.92 25.12 26.33 130.33  2.23 DH
pfout <- polyFit(vert,2,use='glm')
newx <- vert[1,-7]
newx[1] <- 30  # what if V1 were only 30 for case 1?
newx
#   V1    V2    V3    V4    V5    V6
# 1 30 22.55 39.61 40.48 98.67 -0.25
predict(pfout,newx)
# [1] "NO"
```
Forward stepwise regression is also available with `FSR` which also accepts polynomial degree and interaction as inputs. 

```r
out <- FSR(iris)

set seed to -162982340.

The dependent variable is 'Species' which will be treated as multinomial.  The data contains 150 observations (N_train == 113 and N_test == 37), which were split using seed -162982340. The data contains 4 continuous features and 2 dummy variables. Between 6 and 105 models will be estimated. Each model will add a feature, which will be included in subsequent models if it explains at least an additional 0.01 of variance out-of-sample (after adjusting for the additional term on [0, 1]).

Multinomial models will be fit with 'setosa' (the sample mode of the training data) as the reference category.


beginning Forward Stepwise Regression...


# weights:  9 (4 variable)
initial  value 124.143189 
iter  10 value 66.077933
iter  20 value 65.630757
iter  30 value 65.615252
final  value 65.614681 
converged



The added feature WAS accepted into model 1 

(training) AIC: 139.2294 
(training) BIC: 144.6841 
(test) classification accuracy: 0.7027027 

######### (output abbreviated) ##########

The output has a data.frame out$models that contains measures of fit and information about each model, such as the formula call. The output is also a nested list such that if the output is called 'out', out$model1, out$model2, and so on, contain further metadata. The predict method will automatically use the model with the best validated fit but individual models can also be selected like so:

predict(out, newdata = Xnew, model_to_use = 3) 

```
`FSR()` contains a handful of parameters which make the function more or less 'optimistic' about estimating new models. `threshold_include` sets the minimum improvment on the best model to include new features (default 0.01 in adjusted R^2^ for continuous outcomes and accuracy for multinomial outcomes, with the same adjustment applied). `threshold_estimate`, is the treshold to keep adding additional features on the same scale. For categorical outcomes, a linear probability model can also be estimated via Ordinary Least Squares for speed.

```r
out <- FSR(iris, linear_estimation=TRUE)
```

