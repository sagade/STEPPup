\name{OptimalCutPoint}
\alias{OptimalCutPoint}
\title{Cut point search via cross-validation}
\usage{
OptimalCutPoint(response, var, treat, covars = NULL, data, num.cuts = 10000,
  K = 10, repeats = 5, cores = 1, seed = 123)
}
\arguments{
  \item{response}{the name of the binary response variable}

  \item{var}{the name of continuous variable}

  \item{treat}{the name of the treatment variable}

  \item{covars}{optional vector with covariate names}

  \item{data}{the data frame with the data}

  \item{num.cuts}{the number of cuts to be tested. If the
  number of cuts is larger than the number of unique
  values...}

  \item{K}{the number of folds of the cross-validation}

  \item{repeats}{the nu}

  \item{cores}{}

  \item{seed}{the seed to be used for cross-validation}
}
\value{
matrix with cut points, median AUC and AUCs for every
repeat
}
\description{
Searches an optimal cut-point for a continuous variable
'var' in interaction with treatmen 'treat' for a binary
response. The cut point is optimized via cross-validated
AUC
}

