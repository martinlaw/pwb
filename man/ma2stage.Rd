% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ma2stage.R
\name{ma2stage}
\alias{ma2stage}
\title{ma2stage}
\usage{
ma2stage(N = 2428, theta0 = 0.2, theta1 = 0.15, nsims = 1e+05, n.studies = 4)
}
\arguments{
\item{N}{Maximum sample size (equal for all trials)}

\item{theta0}{Response probability for control group}

\item{theta1}{Response probability for experimental group}

\item{nsims}{number of simulations (default 1e5)}

\item{n.studies}{number of non-adaptive studies in meta-analysis}
}
\value{
Data frame of results
}
\description{
Examines bias for a 2-stage binary outcome design WRT a meta-analysis of non-adaptive designs
}
