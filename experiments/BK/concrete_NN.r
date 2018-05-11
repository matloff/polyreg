# Data obtained from: https://archive.ics.uci.edu/ml/datasets/Concrete+Compressive+Strength
# Citation to include in our manuscript: Yeh IC. Modeling of strength of high performance concrete using arti cial neural networks. Cement and Concrete Research. 1998; 28:1797-1808.


## Load in libraries
library(neuralnet)

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

## 75:25 split for training/test
concrete_train <- concrete_normalized[1:773, ]
concrete_test <- concrete_normalized[774:1030, ]

### Train multilayer feedforward neural net with five hidden nodes
concrete_model <- neuralnet(concrete_compressive_strength ~ cement + blast_furnace_slag
                            + fly_ash + water + superplasticizer + coarse_aggregate + fine_aggregate + age,
                            data = concrete_train, hidden = 5)

## Produce a network topology diagram to get a peek inside the black box of our NN with five hidden nodes
plot(concrete_model)
### --- Sum of Squared Errors (SSE): 1.644409 --- ###

## Evaluate model performance on test data
model_results <- compute(concrete_model, concrete_test[1:8])
predicted_concrete_compressive_strength <- model_results$net.result
cor(predicted_concrete_compressive_strength, concrete_test$concrete_compressive_strength)
### --- correlation is 0.6078374271 for five hidden nodes --- ###
