\name{STEPPSurvModel}
\alias{STEPPSurvModel}
\title{STEPPSurvModel}
\usage{
STEPPSurvModel(response, treat, covars = NULL, data, timepoint = 60, alpha)
}
\arguments{
  \item{response}{the response, must be a Surv object}

  \item{treat}{the name of the treatment variable in the
  data}

  \item{covars}{optinal names of covariates}

  \item{data}{the data frame used for evalution}

  \item{the}{timepoint for evaluating the KM estimate, only
  used if no treatment is given}

  \item{alpha}{the alpha level used for confidence
  intervals}
}
\description{
A helper function for the STEPP analysis in case of a
surival response: if a treatment is given a Cox model is
used with optional covariates. If no treatment is given a
survival rate (or more precisely the survival probability)
is estimated using the Kaplan-Meier estimate at a given
timepoint
}

