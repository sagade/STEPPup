\name{ComputeAUC}
\alias{ComputeAUC}
\title{Helper function to compute the AUCs from a prediction matrix}
\usage{
ComputeAUC(predictions, labels, correct = T)
}
\arguments{
  \item{predictions}{matrix with predictions, rows are the
  individuals, cols the repeats}

  \item{labels}{true class labels}

  \item{correct}{}
}
\value{
vector with AUCs, length equal to the number of cols
}
\description{
Helper function to compute the AUCs from a prediction
matrix
}

