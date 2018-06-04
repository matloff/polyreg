# Data obtained from: https://archive.ics.uci.edu/ml/datasets/Concrete+Compressive+Strength
# Citation to include in our manuscript: Yeh IC. Modeling of strength of high performance concrete using arti cial neural networks. Cement and Concrete Research. 1998; 28:1797-1808.


## Load in libraries
library(polyreg)

## Set seed for reproducibility purposes
set.seed(777)

## Load in the already randomly sorted dataset (response variable ("concrete_compressive_strength") is conveniently already in the last column -- in prep for polyreg)
concrete <- read.csv("concreteStrengthData.csv")

## Create a normalization function (since many of the 8 features that predict "concrete_compressive_strength" are all over the place range-wise and are not normally distributed)
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

## Apply the normalization function
concrete_normalized <- as.data.frame(lapply(concrete, normalize))

## polyreg magic
f <- polyFit(concrete_normalized,2,2,"lm",TRUE, 0.9)
pred <- predict(f,concrete_normalized[,-ncol(concrete_normalized)])

## moment of truth...
cor(pred, concrete_normalized$concrete_compressive_strength)
### --- correlation is 0.8693011328, much better than an NN with 5 hidden layers (0.6078374271)--- ###