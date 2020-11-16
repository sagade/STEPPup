\name{STEPP}
\alias{STEPP}
\title{STEPP analysis}
\usage{
STEPP(response, var, treat, covars = NULL, data = parent.frame(),
  mode = c("sw", "to"), alpha = 0.05, adapt.ci = c("none", "normal",
  "bonferroni"), log = F, timepoint, ...)
}
\arguments{
  \item{response}{the response of the model, either a Surv
  object or the name of a factor variable in 'data'}

  \item{var}{the name of the continuous variable}

  \item{treat}{the name of the treatment variable, either a
  numerical vector or a factor. First level is fail, second
  is success}

  \item{covars}{other variables to include in the model.
  Default to NULL}

  \item{data}{the data frame the variables are evaluated
  in}

  \item{mode}{the method to compute subrgoups. Either 'sw'
  for sliding-window or 'to' for tail-oriented}

  \item{alpha}{the alpha for testing. Default to 0.05}

  \item{adapt.ci}{indicates whether a pointwise confidence
  iterval is drawn or an confidence band}

  \item{log}{indicates wheter the Odds-ratios/Hazard ratios
  will be given at a log scale or not. Success rates
  (binary outcome) and survival rates are not affected.
  Default to FALSE}

  \item{...}{further arguments to 'STEPPSliding' or
  'STEPPTail' depending on 'mode'}
}
\value{
object from class 'STEPP'
}
\description{
STEPP analysis
}

