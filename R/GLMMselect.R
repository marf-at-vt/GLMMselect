

bernoulli_pseudolikelihood <- function(Y, X, Beta, Z, Alpha){
  exp_value <- exp(X%*%Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))
  mu <- exp_value/(1+exp_value)
  V <- drop(exp_value/(1+exp_value)^2)
  inv_V_vector <- 1/V
  y_star <- inv_V_vector*(Y-mu)+log(exp_value)
  return(list(mu=mu, V=V, inv_V_vector=inv_V_vector, y_star=y_star))
}

poisson_pseudolikelihood <- function(Y, X, Beta, Z, Alpha, offset){
  if(is.null(offset)){
    offset = 1
  }
  mu <- exp(X%*%Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))*offset
  V <- drop(mu)
  inv_V_vector <- 1/V
  y_star <- inv_V_vector*(Y-mu)+log(mu)-log(offset)
  return(list(mu=mu, V=V, inv_V_vector=inv_V_vector, y_star=y_star))
}

pseudo_likelihood <- function(Y, X, Sigma, Z, family, offset=NULL){
  # pseudo likelihood method, return estimates of parameters, adjusted observations and ...
  # Y: original observations
  # X: covariates
  # Sigma (covariance matrix for random effects): list
  # Z (design matrix for random effects): list

  # family & default link function
  # bernoulli	(link = "logit")
  # poisson	(link = "log")

  # offset: num of measurements for each individual (for poisson distribution)

  X1 <- cbind(1,X)
  if(family == "bernoulli"){
    glmfit <- stats::glm(Y~X, family = "binomial")
  }else{
    glmfit <- stats::glm(Y~X, family = family)
  }

  Beta <- glmfit$coefficients
  n <- length(Y)
  n_rf <- length(Sigma)
  Alpha <- list()
  Kappa <- list()
  Kappa_temp <- list()
  for(i_rf in 1:n_rf){
    Alpha[[i_rf]] <- matrix(rep(0, nrow(Sigma[[i_rf]])), ncol = 1)
    Kappa[[i_rf]] <- 0
    Kappa_temp[[i_rf]] <- 0
  }
  Beta_temp <- matrix(0, nrow = ncol(X1), ncol = 1)

  while(sum(mapply(function(x,y){abs(x-y)>0.0001},Kappa, Kappa_temp))>=1 | max(abs(Beta-Beta_temp))>0.0001){

    Kappa_temp <- Kappa
    Beta_temp <- Beta

    if(family=="bernoulli"){
      par_est <- bernoulli_pseudolikelihood(Y=Y,X=X1, Beta=Beta, Z=Z, Alpha=Alpha)
    }
    if(family=="poisson"){
      par_est <- poisson_pseudolikelihood(Y=Y,X=X1, Beta=Beta, Z=Z, Alpha=Alpha, offset=offset)
    }

    mu <- par_est$mu
    V <- par_est$V
    inv_V_vector <- par_est$inv_V_vector
    y_star <- par_est$y_star

    logl <- function(x){
      # log posterior of kappa
      H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, x, Sigma, Z)), ncol=n, nrow=n)
      diag(H) <- diag(H)+inv_V_vector
      eign <- eigen(H, symmetric = TRUE)
      P <- eign$vectors
      D <- eign$values
      D_inv <- 1/D
      D_inv[!is.finite(D_inv)] <- 0
      X1.tilde <- t(P)%*%X1

      re.tilde <- t(P)%*%y_star-X1.tilde%*%Beta

      l <- (+0.5*sum(log(D_inv))
            -0.5*determinant(t(X1.tilde)%*%(D_inv*X1.tilde), logarithm=TRUE)$modulus[1]
            -0.5*sum(re.tilde^2*D_inv))
      return(-l)
    }

    ts <- stats::optim(par=rep(0.1,n_rf), fn=logl, lower = rep(0,n_rf), upper = rep(50,n_rf), method = "L-BFGS-B",hessian = T)
    Kappa <- as.list(ts$par)
    Hessian <- ts$hessian

    H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, Kappa, Sigma, Z)), ncol=n, nrow=n)
    diag(H) <- diag(H)+inv_V_vector
    eign <- eigen(H, symmetric = TRUE)
    D <- eign$values
    P <- eign$vectors
    D_inv <- 1/D
    D_inv[!is.finite(D_inv)] <- 0

    X1.tilde <- t(P)%*%X1
    y.tilde <- t(P)%*%y_star
    #update beta and alpha
    xdx <- t(X1.tilde)%*%(D_inv*X1.tilde)
    eign <- eigen(xdx, symmetric = TRUE)
    D_xdx <- eign$values
    P_xdx <- eign$vectors
    D_inv_xdx <- 1/D_xdx
    D_inv_xdx[!is.finite(D_inv_xdx)] <- 0
    inv_xdx <- P_xdx%*%(D_inv_xdx*t(P_xdx))
    Beta <- inv_xdx%*%t(X1.tilde)%*%matrix(D_inv*y.tilde,ncol = 1)
    re.tilde <- y.tilde-X1.tilde%*%Beta

    Alpha <- mapply(function(x,y,z){x*y%*%t(z)%*%P%*%matrix(D_inv*re.tilde,ncol=1)}, Kappa, Sigma, Z, SIMPLIFY = FALSE)

  }


  if(family=="bernoulli"){
    par_est <- bernoulli_pseudolikelihood(Y=Y,X=X1, Beta=Beta, Z=Z, Alpha=Alpha)
  }
  if(family=="poisson"){
    par_est <- poisson_pseudolikelihood(Y=Y,X=X1, Beta=Beta, Z=Z, Alpha=Alpha, offset=offset)
  }

  mu <- par_est$mu
  V <- par_est$V
  inv_V_vector <- par_est$inv_V_vector
  y_star <- par_est$y_star

  est <- c(as.vector(Beta)[-1], unlist(Kappa))
  sd <- c(sqrt(diag(inv_xdx))[-1], sqrt(1/diag(Hessian)))

  return(list(y_star=y_star, inv_v=inv_V_vector, est=est, sd=sd))
}



FBF <- function(y_star, Xc, inv_v, Zc, Sigmac, prior,b, K){
  # calculate posterior prob by fractional Bayes factor method, return log of posterior prob
  #y_star: adjusted observations
  #Xc: covariates in the model c
  #inv_v: inverse matrix of V
  # Zc: design matrix for random effects in the model c, list
  # Sigmac: covariance matrix for random effects in the model c, list
  #prior: "AR", "HC" prior
  #b: train size
  #K: num of candidate covariates

  kc <- ncol(Xc)                 #num of covariates in model c
  Xc <- cbind(1, Xc)             #X matrix in model c
  N <- nrow(y_star)
  q <- length(Sigmac)

  logmarginallikelihood <- function(a,b,prior){
    #log marginal likelihood of variance component, tau
    #a:list of transform of tau, b: train size

    H_t <- matrix(rowSums(mapply(function(x,y,z){exp(x)*z%*%y%*%t(z)}, a, Sigmac, Zc)), ncol=N, nrow=N)
    diag(H_t) <- diag(H_t)+inv_v

    eign <- eigen(H_t,symmetric = TRUE)
    eign.value <- eign$values
    eign.value[abs(eign$values)<10^(-10)] <- 0
    P <- eign$vectors
    D <- eign.value^(-1)
    D[!is.finite(D)] <- 0
    inv_H_t <- P%*%(D*t(P))  # inverse of H_t

    XHX <- t(Xc)%*%inv_H_t%*%Xc
    if (kc==0){
      inv_XHX <- 1/XHX
    } else{
      eign <- eigen(XHX,symmetric = TRUE)
      eign.value <- eign$values
      eign.value[abs(eign$values)<10^(-10)] <- 0
      P <- eign$vectors
      D <- eign.value^(-1)
      D[!is.finite(D)] <- 0
      inv_XHX <- P%*%(D*t(P))  # inverse of XHX

    }
    if(prior == "HC"){
      r <- -(b/2*determinant(inv_H_t)$modulus[1]+determinant(inv_XHX)$modulus[1]/2-(kc+1)/2*log(b)
             +sum(sapply(a, function(x){x-log(exp(x)+1)}))
             +b/2*drop(t(y_star)%*%(inv_H_t%*%Xc%*%inv_XHX%*%t(Xc)%*%inv_H_t-inv_H_t)%*%y_star))
    }

    if(prior == "AR"){
      r <- -(b/2*determinant(inv_H_t)$modulus[1]+determinant(inv_XHX)$modulus[1]/2-(kc+1)/2*log(b)
             +sum(sapply(a, function(x){x-2*log(exp(x)/2+1)}))
             +b/2*drop(t(y_star)%*%(inv_H_t%*%Xc%*%inv_XHX%*%t(Xc)%*%inv_H_t-inv_H_t)%*%y_star))
    }

    return(r)
  }


  logmarginallikelihood_null <- function(b){
    #log marginal likelihood of variance component, tau
    # b: train size

    V <- diag(1/inv_v,nrow = N, ncol = N)
    XHX <- t(Xc)%*%V%*%Xc
    if (kc==0){
      inv_XHX <- 1/XHX
    } else{
      eign <- eigen(XHX,symmetric = TRUE)
      eign.value <- eign$values
      eign.value[abs(eign$values)<10^(-10)] <- 0
      P <- eign$vectors
      D <- eign.value^(-1)
      D[!is.finite(D)] <- 0
      inv_XHX <- P%*%(D*t(P))  # inverse of XHX

    }

    r <- (b/2*determinant(V)$modulus[1]+determinant(inv_XHX)$modulus[1]/2-(kc+1)/2*log(b)
          +b/2*drop(t(y_star)%*%(V%*%Xc%*%inv_XHX%*%t(Xc)%*%V-V)%*%y_star))

    return(r)
  }

  laplaceapprox <- function(b,prior){
    # use laplace approx to get integration of tau
    par <- rep(-2,q)
    lower <- rep(-10,q)
    upper <- rep(5,q)
    op <- stats::optim(par=par, fn=logmarginallikelihood, b=b, prior=prior, lower = lower, upper = upper, method = "L-BFGS-B", hessian = TRUE)
    ap <- -op$value-q/2*log(2*pi)-log(det(op$hessian))/2
    return(ap)
  }

  if(q != 0){
    qc_N <- laplaceapprox(b=1, prior=prior)    # log numerator of fractional integrated likelihood
    qc_D <- laplaceapprox(b=b, prior=prior)    # log denominator of fractional integrated likelihood
  }else{
    qc_N <- logmarginallikelihood_null(1)
    qc_D <- logmarginallikelihood_null(b)
  }

  prior_Mc <- 1/choose(K, kc)  #proportion of prior prob
  log_posterior_Mc <- qc_N-qc_D+log(prior_Mc) #log of posterior of model c

  return(log_posterior_Mc)

}

#' GLMMselect: Bayesian model selection method for generalized linear mixed models

#' @importFrom stats optim

#' @param Y A numeric vector of observations.
#' @param X A matrix of covariates.
#' @param Sigma A list of covariance matrices for random effects.
#' @param Z A list of design matrices for random effects.
#' @param family A description of the error distribution to be used in the model.
#' @param prior The prior distribution for variance component of random effects.
#' @param offset This can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of observations.
#' @param NumofModel The number of models with the largest posterior probabilities being printed out.
#' @param pip_fixed The cutoff that if the posterior inclusion probability of fixed effects is larger than it, the fixed effects will be included in the best model.
#' @param pip_random The cutoff that if the posterior inclusion probability of random effects is larger than it, the random effects will be included in the best model.
#' @returns A list of the indices of covariates and random effects which are in the best model.
#' @examples
#'
#' library(GLMMselect)
#' \donttest{
#' data("Y");data("X");data("Z");data("Sigma")
#' Model_selection_output <- GLMMselect(Y=Y, X=X, Sigma=Sigma,
#'                          Z=Z, family="poisson", prior="AR", offset=NULL)
#' }
#' @export
GLMMselect <- function(Y, X, Sigma, Z, family, prior, offset=NULL, NumofModel=10, pip_fixed=0.5, pip_random=0.5){

  family <- tolower(family)
  if(sum(family %in% c("poisson","bernoulli")) == 0){
    stop("Family must be either Poisson or Bernoulli.")
  }


  if(!is.numeric(Y)){
    stop("Y has to be numeric.")
  }

  if(!is.matrix(X)){
    stop("X has to be a matrix object.")
  }

  if(!is.numeric(X)){
    stop("X has to contain numeric values.")
  }

  if(sum(prior %in% c("AR","HC")) == 0){
    stop("Prior must be either AR or HC.")
  }

  if(family == "bernoulli"){
    if(!is.null(offset)){
      stop("Offset is only for Poisson data.")
    }
  }

  if(!is.list(Sigma)){
    stop("Sigma has to be a list object.")
  }

  if(!is.list(Z)){
    stop("Z has to be a list object.")
  }

  if(sum(mapply(is.numeric, Sigma))==0){
    stop("Sigma must be numeric.")
  }

  if(sum(mapply(is.numeric, Z))==0){
    stop("Z must be numeric.")
  }

  if(family == "bernoulli"){
    family = "binominal"
  }

  N <- length(Y)
  X1 <- cbind(1,X)
  K <- ncol(X)
  p <- K+1
  b <- (p+1)/N #FBF train size
  Q <- length(Sigma)
  PL_est <- pseudo_likelihood(Y=Y, X=X, Sigma=Sigma, Z=Z, family=family, offset=offset)

  postprob <- matrix(NA, nrow = 2^Q, ncol = 2^K)
  PosteriorProb <- matrix(NA, nrow = 2^Q*2^K, ncol = (K+Q+1))
  colnames(PosteriorProb) <- c(paste("x",1:K,sep = ""),paste("r",1:Q,sep = ""),"MPP")


  binary_fixed <- rep(list(0:1), K)
  binary_fixed <- as.matrix(expand.grid(binary_fixed))

  binary_random <- rep(list(0:1), Q)
  binary_random <- as.matrix(expand.grid(binary_random))

  for(q in 1:2^Q){
    for(k in 1:2^K){
      Xc <- X[, which(binary_fixed[k,]==1),drop=FALSE]
      Sigmac <- Sigma[which(binary_random[q,]==1)]
      Zc <- Z[which(binary_random[q,]==1)]
      postprob[q,k] <- FBF(y_star=PL_est$y_star, Xc=Xc, inv_v=PL_est$inv_v, Zc=Zc, Sigmac=Sigmac, prior=prior, b=b, K=K)
    }
  }

  maxx <- max(postprob)
  postprob <- exp(postprob-maxx)/sum(exp(postprob-maxx))

  j=1
  for(q in 1:2^Q){
    for(k in 1:2^K){
      PosteriorProb[j,1:K] <- binary_fixed[k,]
      PosteriorProb[j,(K+1):(K+Q)] <- binary_random[q,]
      PosteriorProb[j,K+Q+1] <- postprob[q,k]
      j = j+1
    }
  }


  ###output###

  #table of models' posterior probabilities
  PosteriorProb <- PosteriorProb[order(PosteriorProb[,K+Q+1],decreasing=TRUE),]
  PosteriorProb <- as.data.frame(PosteriorProb)
  PosteriorProb$MPP <- round(PosteriorProb$MPP,3)

  indices <- apply(PosteriorProb[,1:(K+Q)],2,function(x){which(x==1)})
  margin_prob <- apply(indices,2,function(x){sum(PosteriorProb[x,K+Q+1])})
  #table for fixed effects
  FixedEffect <- matrix(NA, ncol = 3, nrow = K)
  rownames(FixedEffect) <- paste("x",1:K,sep = "")
  colnames(FixedEffect) <- c("Est","SD","PIP")
  FixedEffect[,1] <- PL_est$est[1:K]
  FixedEffect[,2] <- PL_est$sd[1:K]
  FixedEffect[,3] <- margin_prob[1:K]
  FixedEffect[,1:2] <- round(FixedEffect[,1:2],3)
  FixedEffect[,3] <- round(FixedEffect[,3],3)

  #table for random effects
  RandomEffect <- matrix(NA, ncol = 3, nrow = Q)
  rownames(RandomEffect) <- paste("r",1:Q,sep = "")
  colnames(RandomEffect) <- c("Est","SD","PIP")
  RandomEffect[,1] <- PL_est$est[(1+K):(Q+K)]
  RandomEffect[,2] <- PL_est$sd[(1+K):(Q+K)]
  RandomEffect[,3] <- margin_prob[(1+K):(Q+K)]
  RandomEffect[,1:2] <- round(RandomEffect[,1:2],3)
  RandomEffect[,3] <- round(RandomEffect[,3],3)


  if(NumofModel <= (2^K*2^Q)){
    PosteriorProb <- PosteriorProb[1:NumofModel,]
  }


  #bestmodel
  BestModel <- list()
  fix_position <- which(margin_prob[1:K]>pip_fixed)
  random_position <- which(margin_prob[(1+K):(Q+K)]>pip_random)
  names(fix_position)<-NULL
  names(random_position)<-NULL
  BestModel$covariate_inclusion <- fix_position
  BestModel$random_effect_inclusion <- random_position


  return(list(BestModel=BestModel, PosteriorProb=PosteriorProb, FixedEffect=FixedEffect, RandomEffect=RandomEffect))

}




