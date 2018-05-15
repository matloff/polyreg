## Load in the already randomly sorted dataset (response variable ("concrete_compressive_strength") is conveniently already in the last column -- in prep for polyreg)
concrete <- read.csv("experiments/BK/concreteStrengthData.csv")

library(kerasformula)
out <- kms(concrete_compressive_strength ~., concrete, pTraining=0.9, validation_split = 0, seed=777)
out$R2_predictions

# 0.546152
# implies correlation between out of sample prediction and y_test is 0.739021
# between NN and polyreg at default kms settings...
