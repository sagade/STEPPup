\name{STEPPTailDown}
\alias{STEPPTailDown}
\title{Helper function to get a tail-orinted (downwards: starting from all to few patients) subpopulation pattern}
\usage{
STEPPTailDown(x, g)
}
\arguments{
  \item{x}{vector with covariate}

  \item{g}{desired number of groups}
}
\value{
list with k elements, each holding the indices of the
elements from x belonging to that group
}
\description{
Helper function to get a tail-orinted (downwards: starting
from all to few patients) subpopulation pattern
}
\examples{
set.seed(123)
x <- runif(200)
STEPPTailDown(x, 10)
}

