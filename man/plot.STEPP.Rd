\name{plot.STEPP}
\alias{plot.STEPP}
\title{Plot function for an STEPP object}
\usage{
\method{plot}{STEPP}(x, type = c("coef", "group", "coef_group", "coef_box",
  "coef_density"), coef.xaxis = c("medians", "groups"), coef.ylab, plot = T,
  ...)
}
\arguments{
  \item{x}{STEPP object as returned by STEPP function}

  \item{type}{the type of the plot to be produced: can be
  'coef' for just the coefficients, 'group' for just the
  groupings, 'coef_group' for both plots aligned, or
  'coef_box' and 'coef_density' for the coefficient plot
  with either a boxplot or density plot displaying the
  distribution of the variable of interest}

  \item{coef.xaxis}{type of the x-axis of the coefficient
  plot, either 'medians' for a continuous scale using the
  medians of the variable of interest in each subgroup.
  'groups' uses a discrete (and equidistant) x-axis with
  group names (limits of the variable of interest in that
  subgroup). The parameter has only an effect when type is
  'coef', for combined plots it is set automatically.}

  \item{coef.ylab}{the y-axis label for the coefficient
  plot. It missing the label will be guessed based on the
  type of the analysis: cox, logistic, survrate or rate}

  \item{plot}{should the resulting plot be plotted or only
  returned ?}

  \item{...}{further arguments to STEPPCoefPlot or
  STEPPGroupPlot if type is 'group'}
}
\description{
Plot function for an STEPP object
}

