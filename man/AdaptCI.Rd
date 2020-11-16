\name{AdaptCI}
\alias{AdaptCI}
\title{Bonferroni adaption of the confidence interval}
\usage{
AdaptCI(coefs, alpha, k = nrow(coefs))
}
\arguments{
  \item{coefs}{matrix with coefficients and confidence
  intervals. The first column holds the coeffecients, the
  second the lower and the third the upper bound of the
  interval}

  \item{alpha}{the \alpha from the CI}

  \item{k}{number of subgroups. Default to number of rows
  from the coefs matrix}
}
\description{
Bonferroni adaption of the confidence interval
}

