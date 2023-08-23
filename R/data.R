#' The vector of observations
#'
#' The data was simulated under a Poisson Generalized linear mixed model with four candidate covariates, a vector of spatial random effects, and a vector of overdispersion random effects.
#' The matrix of candidate covariates, the design matrices for two types of random effects, and the covariance matrices for two types of random effects are provide in the package.
#'
#' @format ## `Y`
#' A vector with 100 observations.
"Y"

#' The matrix of candidate covariates
#'
#' This is a matrix with four candidate covariates. The covariates in the first two columns are in the true model.
#'
#' @format ## `X`
#' A matrix with 100 observations and 4 covariates.
"X"

#' A list of design matrices
#'
#' This is a list of two design matrices for two types of random effects.
#' The first component is for spatial random effects. The second componet is for overdispersion random effects.
#'
#' @format ## `Z`
#' A list of design matrices for random effects.
"Z"

#' A list of covariance matrices
#'
#' This is a list of two covariance matrices for two types of random effects.
#' The first component is for spatial random effects. The second componet is for overdispersion random effects.
#'
#' @format ## `Sigma`
#' A list of covariance matrices for random effects.
"Sigma"







#' The vector of observations for lip cancer case study
#'
#' The observations about the observed number of lip cancer cases in each Scotland county between 1975 and 1980.
#'
#' @format ## `lipcancer_Y`
#' A vector with 56 observations.
"lipcancer_Y"

#' The vector of priori information for lip cancer case study
#'
#' The expected number of lip cancer cases in each Scotland county between 1975 and 1980.
#'
#' @format ## `lipcancer_offset`
#' A vector
"lipcancer_offset"


#' The matrix of candidate covariates
#'
#' This is a matrix with one candidate covariates, which is the proportion of population engaged in agriculture, fishing, or forestry.
#'
#' @format ## `lipcancer_X`
#' A matrix with 56 observations and 1 covariate.
"lipcancer_X"

#' A list of design matrices
#'
#' This is a list of two design matrices for two types of random effects.
#' The first component is for spatial random effects. The second componet is for overdispersion random effects.
#'
#' @format ## `lipcancer_Z`
#' A list of design matrices for random effects.
"lipcancer_Z"

#' A list of covariance matrices
#'
#' This is a list of two covariance matrices for two types of random effects.
#' The first component is for spatial random effects. The second componet is for overdispersion random effects.
#'
#' @format ## `lipcancer_Sigma`
#' A list of covariance matrices for random effects.
"lipcancer_Sigma"

