\name{asilb}
\alias{asilb}
\alias{asilb.formula}
\alias{asilb.default}
\alias{asilb.data.frame}
\alias{asilb.matrix}
\title{(Adjacent-categories) Simple Isotonic LogitBoost}
\description{Train and predict logitboost-based classification algorithm using isotonic regression (decision stumps for no monotone features) as weak learners, based on the adjacent-categories logistic model (see Agresti (2010)). For full details on this algorithm, see Conde et al. (2020).}
\usage{
asilb(xlearn, \dots)

\method{asilb}{formula}(formula, data, \dots)

\method{asilb}{default}(xlearn, ylearn, xtest = xlearn, mfinal = 100, 
monotone_constraints = rep(0, dim(xlearn)[2]), prior = NULL, \dots)
}
\arguments{
  \item{formula}{ A formula of the form \code{groups ~ x1 + x2 + ...}. That is, the response is the class variable and the right hand side specifies the explanatory variables. }
  \item{data}{ Data frame from which variables specified in \code{formula} are to be taken. }
  \item{xlearn}{ (Required if no formula is given as the principal argument.) A data frame or matrix containing the explanatory variables. }
  \item{ylearn}{ (Required if no formula is given as the principal argument.) A numeric vector or factor with numeric levels specifying the class for each observation. }
  \item{xtest}{ A data frame or matrix of cases to be classified, containing the features used in \code{formula} or \code{xlearn}. }
  \item{mfinal}{ Number of iterations of the algorithm. }
  \item{monotone_constraints}{ Numerical vector consisting of 1, 0 and -1, its length equals the number of features in \code{xlearn}. 1 is increasing, -1 is decreasing and 0 is no constraint. }
  \item{prior}{ The prior probabilities of class membership. If unspecified, equal prior probabilities are used. If present, the probabilities must be specified in the order of the factor levels. }
  \item{\dots}{ Arguments passed to or from other methods. }
}
\value{
A list containing the following components:
  \item{call}{ The (matched) function call. }
  \item{trainset}{ Matrix with the training set used (first columns) and the class for each observation (last column). }
  \item{prior}{ Prior probabilities of class membership used. }
  \item{apparent}{ Apparent error rate. }
  \item{mfinal}{ Number of iterations of the algorithm. }
  \item{loglikelihood}{ Log-likelihood. }
  \item{posterior}{ Posterior probabilities of class membership for \code{xtest} set. }
  \item{class}{ Labels of the class with maximal probability for \code{xtest} set. }
}
\references{
Agresti, A. (2010). Analysis of Ordinal Categorical Data, 2nd edition. John Wiley and Sons. New Jersey.

Conde, D., Fernandez, M. A., Rueda, C., and Salvador, B. (2020). Isotonic boosting classification rules. \emph{Advances in Data Analysis and Classification}, 1-25.
}
\author{
David Conde
}
\note{
This function may be called using either a formula and data frame, or a data frame and grouping variable, or a matrix and grouping variable as the first two arguments. All other arguments are optional.

Classes must be identified, either in a column of \code{data} or in the \code{ylearn} vector, by natural numbers varying from 1 to the number of classes. The number of classes must be greater than 1.

If there are missing values in either \code{data}, \code{xlearn} or \code{ylearn}, corresponding observations will be deleted.
}
\seealso{
\code{\link{amilb}}, \code{\link{csilb}}, \code{\link{cmilb}}
}
\examples{
data(motors)
table(motors$condition)
##  1  2  3  4 
## 83 67 70 60 

## Let us consider the first three variables as predictors
data <- motors[, 1:3]
grouping = motors$condition
## 
## Lower values of the amplitudes are expected to be 
## related to higher levels of damage severity, so 
## we can consider the following monotone constraints
monotone_constraints = rep(-1, 3)

set.seed(7964)
values <- runif(dim(data)[1])
trainsubset <- values < 0.2
obj <- asilb(data[trainsubset, ], grouping[trainsubset], 
               data[-trainsubset, ], 50, monotone_constraints)

## Apparent error
obj$apparent
## 4.761905

## Error rate
100*mean(obj$class != grouping[-trainsubset])
## 14.69534
}
\keyword{ multivariate }
\keyword{ classif }