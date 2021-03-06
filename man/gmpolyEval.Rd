% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval.R
\name{gmpolyEval}
\alias{gmpolyEval}
\title{Evaluation of a multivariate polynomial}
\usage{
gmpolyEval(pol, x)
}
\arguments{
\item{pol}{a \code{\link{gmpoly}} object}

\item{x}{either a \code{\link[gmp]{bigq}} vector, or a 
\code{\link[gmp]{bigq}} matrix; in the later case, the evaluation is 
performed for each row of the matrix}
}
\value{
A \code{\link[gmp]{bigq}} number or vector.
}
\description{
Evaluates a \code{gmpoly} multivariate polynomial for given 
  values of the variables.
}
\examples{
library(gmpoly)
library(gmp)
pol <- gmpoly("5/2 x^(2,3) + 3 x^(1,1)")
gmpolyEval(pol, as.bigq(c(1, 1)))
x <- rbind(
  t(as.bigq(c(1, 1))),
  t(as.bigq(c(3, 4), c(4, 3)))
)
gmpolyEval(pol, x)
}
