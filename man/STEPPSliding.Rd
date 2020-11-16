\name{STEPPSliding}
\alias{STEPPSliding}
\title{Helper function to get a sliding-window subpopulation pattern}
\usage{
STEPPSliding(x, n2, n1)
}
\arguments{
  \item{x}{vector with covariate}

  \item{n2}{size of each window (number of patients in one
  window)}

  \item{n1}{min difference between windows ('slide width')}
}
\value{
list with k elements, each holding the indices of the
elements from x belonging to that group
}
\description{
Helper function to get a sliding-window subpopulation
pattern
}
\examples{
set.seed(123)
x <- runif(200)
STEPPSliding(x, n2=10, n1=5)
}
\references{
Bonetti2000_StatInMed
}

