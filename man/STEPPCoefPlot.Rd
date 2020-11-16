\name{STEPPCoefPlot}
\alias{STEPPCoefPlot}
\title{Coef plot from a STEPP analysis}
\usage{
STEPPCoefPlot(coefs, groups = NULL, line = c("smooth", "simple", "none"),
  hline, xlab = "x", ylab = "OR", main = "", log, log.scale = T,
  ylim = NULL)
}
\arguments{
  \item{coefs}{matrix with coefficients and CI}

  \item{line}{for 'smooth' a lowess curve is plotted
  through the coefficients, 'simple' connects all the
  points with a line, 'none' omits the line}

  \item{hline}{intercept for help lines, if NULL no lines
  will be drawn. If missing the log parameter is used to
  draw a line either at 0 (log=TRUE) or 1 (log=FALSE)}

  \item{xlab}{the x axis label}

  \item{xlab}{the y axis label}

  \item{log}{indicates wheter the coefficients are on a log
  scale.coef}

  \item{log.scale}{indicates wheter the y axis should use a
  log scale. Will be ignored if the data are already on a
  log scale}
}
\value{
ggplot object with the coefficient plot
}
\description{
Coef plot from a STEPP analysis
}

