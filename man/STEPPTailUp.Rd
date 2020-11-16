\name{STEPPTailUp}
\alias{STEPPTailUp}
\title{Helper function to get a tail-orinted (upwards: starting from few to all patients) subpopulation pattern}
\usage{
STEPPTailUp(x, g)
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
Helper function to get a tail-orinted (upwards: starting
from few to all patients) subpopulation pattern
}
\examples{
set.seed(123)
x <- runif(200)
STEPPTailUp(x, 10)
}

