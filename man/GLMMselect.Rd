% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GLMMselect.R
\name{GLMMselect}
\alias{GLMMselect}
\title{Bayesian model selection method for generalized linear mixed models}
\usage{
GLMMselect(Y, X, Sigma, Z, family, prior, offset = NULL, NumofModel = 10)
}
\arguments{
\item{Y}{A numeric vector of observations.}

\item{X}{A matrix of covariates.}

\item{Sigma}{A list of covariance matrices for random effects.}

\item{Z}{A list of design matrices for random effects.}

\item{family}{A description of the error distribution to be used in the model.}

\item{prior}{The prior distribution for variance component of random effects.}

\item{offset}{This can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of observations.}

\item{NumofModel}{The number of models with the largest posterior probabilities being printed out.}
}
\value{
A list of the indices of covariates and random effects which are in the best model.
}
\description{
Bayesian model selection method for generalized linear mixed models
}
