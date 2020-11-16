\name{STEPPGroupPlot}
\alias{STEPPGroupPlot}
\title{Helper function to produce a plot that displayes the subpopulations (groups) and number of patients}
\usage{
STEPPGroupPlot(groups, intervals = NULL, xlab = "Subpopulations [No]",
  ylab = "y", main = NULL)
}
\arguments{
  \item{groups}{list with groups}

  \item{intervals}{two-dimensional matrix with lower and
  upper bound of the groups. If none is given (default) the
  names of the group list are parsed}

  \item{xlab}{the x axis label}

  \item{ylab}{the y axis label}

  \item{main}{the plot title. Default to NULL (no title
  set)}
}
\value{
ggplot object with the group plot
}
\description{
Helper function to produce a plot that displayes the
subpopulations (groups) and number of patients
}

