\name{STEPPLogisticModel}
\alias{STEPPLogisticModel}
\title{STEPPLogisticModel}
\usage{
STEPPLogisticModel(response, treat, covars = NULL, data, alpha, verbose = T)
}
\arguments{
  \item{response}{the name of the response in the data}

  \item{treat}{the name of the treatment variable in the
  data}

  \item{covars}{optinal names of covariates}

  \item{data}{the data frame used for evalution}

  \item{alpha}{the alpha level used for confidence
  intervals}
}
\description{
A helper function for the STEPP analysis in case of a
binary response: fits a logistic regression model if a
treatment and optional covars are given. If no treatment is
given a rate is calculated. In the latter case confidence
intervals are calculates based on the biomial distribution.
}

