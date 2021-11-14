# polyreg, an Alternative to Machine Learning Methods

A package to automate formation and evaluation of
multivariate polynomial regression models, especially as an alternative
to neural networks and other machine learning algorithms.  

An important feature is that dummy variables are handled properly, so
that for instance powers of a dummy variable do not exist as duplicates
of the original.

**Note:**  This library is also used in the **qeML** package, with a
convenient, consistent interface, and with extensions such as ridge
polynomial regression.

##Motivation##

In [Polynomial Regression As an Alternative to Neural
Nets](https://arxiv.org/abs/1806.06850), by Cheng, Khomtchouk, Matloff
and Mohanty, 2018, it is shown that dense, feedforward neural networks
are essentially polynomial regression models.  This was extended in
[Towards a Mathematical Framework to Inform Neural Network Modelling via Polynomial Regression](https://www.meta.org/papers/towards-a-mathematical-framework-to-inform-neural/33984736). 
by Morala, Cifuentes, Lillo, and Iñaki Ucar.  The point is then, why go
through the problems of neural networks--convergence, local minima and
so on--when can can work more simply with polynomials?

Of course, it is not quite that simple.  If we start with p variables in
our model, the d-degree polynomial version will have O(p<sup>d</sup>)
variables, which can easily become computationally challenging.
Nevertheless, our experiments have had quite encouraging results.

##Usage##

Other than the various cross-validation functions, the main functions
are **polyfit()** and **predict.polyFit()**.  One can fit either
regression or classification models, with an option to perform PCA for
dimension reduction on the predictors/features.

## Example

Programmer/engineer 2000 Census data, Silicon Valley, built-in to the
package.  


``` r
data(pef)

# model wage income, fitting a degree-2 model
pfout <- polyFit(pef[,c(1,2,3,4,6,5)],2,use='lm')

# predict wage of person like that in row 1, but age 40 and female
newx <- pef[1,-5]
newx$age <- 48
newx$sex <- 1
predict(pfout,newx)  # about $84,330

```

