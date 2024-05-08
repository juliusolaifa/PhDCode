glmm <- function(formula, data, family, computepis=T) {
  modObj <- try(glmmTMB::glmmTMB(formula=formula, data=data, family=family), silent = T)
    
  warning_code <- as.numeric(gsub("[^0-9]", "", modObj$fit$message))
  if(warning_code == 7 || modObj$sdr$pdHess == F) {
    return(NULL)
  }
  
  # Check if modObj is an error
  else if(class(modObj) == "try-error") {
    message("Model fitting failed")
    return(NULL)
  }
  
  else{
    beta <- unname(fixef(modObj)$cond)
    Sigma <- unname(unlist(VarCorr(modObj)$cond))
    phi <- glmmTMB::sigma(modObj)
    p <- unname(glmmTMB::family_params(modObj))
    family <- modObj$modelInfo$family
    vcov  <-  stats::vcov(modObj, full = TRUE)
    logLik <- as.numeric(logLik(modObj))
    group_var <- modObj$modelInfo$grpVar
    nfixed <- length(c(beta,phi))
    nvarcomp_bounded <- length(attr(VarCorr(modObj)$cond[[group_var]], "stddev"))
  
    if(nvarcomp_bounded == 2 && computepis == T) {
      indx <- nfixed+1
      sig12cov <- matrix(c(vcov[indx,indx], vcov[indx,indx+2], 
                           vcov[indx+2,indx], vcov[indx+2,indx+2]),
                         nrow=2,ncol=2)
      pis <- unlist(bivnorm(sig12cov))
    } else{
      pis <- c(0.5, 0.5)
    }
    
    n <- modObj$modelInfo$nobs
  
    result <- list("beta"=beta, "Sigma"=Sigma, "phi"=phi, "p"=p,
                   "family"=family$family, "link"=family$link,
                   "invlink"=family$linkinv,"nfixed" = nfixed,
                   "vcov" = vcov, "Itheta" = n*vcov,
                   "nvarcomp_bounded"=nvarcomp_bounded,
                   "n" = n,  "logLikelihood"=logLik,
                   "modObj"=modObj, "pis"=pis)
    class(result) <- "glmmMod"
    return(result)
  }
}


vpc <- function(fitObj, x=NULL) {
  if(is.null(fitObj)) {
      return(NULL)
  }
  
  if(length(fitObj$beta) == 2 && is.null(x)) {
      stop("x can not be NULL if length(beta) = 2")
  } else if(length(fitObj$beta) == 1) {
      x = 0
  }
  
  # Extract parameters from fitObj
  b0 <- fitObj$beta[1]
  b1 <- if (fitObj$nfixed == 3) fitObj$beta[2] else 0
  phi <- fitObj$phi
  p <- fitObj$p
  n <- fitObj$n
  sig11 <- fitObj$Sigma[1]
  sig12 <- ifelse(fitObj$nvarcomp_bounded == 2, fitObj$Sigma[2], 0)
  sig22 <- ifelse(fitObj$nvarcomp_bounded == 2, fitObj$Sigma[4], 0)
  # Define the vpc function
  vpc_func = get(paste0("vpc_",fitObj$family))
  # Compute the value of vpc
  vpc_value <- switch(fitObj$family,
                      "nbinom2" = vpc_func(b0,b1,phi,sig11,sig12,sig22,x),
                      "tweedie" = vpc_func(b0,b1,phi,sig11,sig12,sig22,p,x)
  )
  
  # Wrapper function for gradient calculation
  vpc_wrapper <- function(params) {
    switch(fitObj$family,
           "nbinom2" = vpc_func(params[1], params[2], params[3], params[4], 
                                params[5], params[6], x),
           "tweedie" = vpc_func(params[1], params[2], params[3], params[4], 
                                params[5], params[6], params[7], x)
    )
  }
  # Parameters vector for gradient calculation
  params <-switch(fitObj$family,
                  "nbinom2" = c(b0, b1, phi, sig11, sig12, sig22),
                  "tweedie" = c(b0, b1, phi, sig11, sig12, sig22, p),
  )
  # Calculate gradient
  gradient <- numDeriv::grad(func = vpc_wrapper, x = params)
  #gradient <- gradient[gradient != 0]
  # Standard Error
  #stderr <- sqrt((gradient %*% fitObj$vcov  %*% gradient)/n)
  # Return both vpc value and gradient
  result <- list("vpc" = vpc_value, "gradient" = gradient, "modObj" = fitObj$modObj,
                 "pis" = fitObj$pis, "vcov" = fitObj$vcov, "Itheta"= fitObj$Itheta, 
                 n =fitObj$n)
  class(result) <- "vpcObj"
  return(result)
}

vpc_nbinom2 <- function(b0, b1, phi, sig11, sig12, sig22, x) {
  mu <- b0 + b1 * x
  Sigma <- matrix(c(sig11, sig12, sig12, sig22), 2, 2)
  Z <- c(1, x)
  sig <- Z %*% Sigma %*% Z
  cmp1 <- (exp(sig) - 1) * exp(2 * mu + sig) # varLogNormal
  cmp2 <- exp(mu + sig / 2)  # meanLogNormal
  cmp3 <- (1 / phi) * exp(2 * mu + 2 * sig)  # 1/phi X 2ndMomentLogNormal
  #cmp3 <- phi * exp(2 * mu + 2 * sig)  # phi X 2ndMomentLogNormal
  return(cmp1 / (cmp1 + cmp2 + cmp3))  # varfunc = mu +mu^2/phi
}

vpc_tweedie <-  function(b0, b1, phi, sig11, sig12, sig22, p, x) {
  mu <- b0 + b1 * x
  Sigma <- matrix(c(sig11, sig12, sig12, sig22), 2, 2)
  Z <- c(1, x)
  sig <- Z %*% Sigma %*% Z
  cmp1 <- (exp(sig) - 1) * exp(2 * mu + sig) # varLogNormal
  cmp2 <- phi* exp(p*mu + 0.5*p^2*sig) #phi X pthMomentLogNormal
  return(cmp1/(cmp1 + cmp2))  # varfunc = phi*mu^p
}


##### Estimating pi's ######
bivnorm <- function(sigma) {
  #print(eigen(sigma)$values)
  bivarnormpdf <- function(x1, x2, sig11, sig22, rho) {
    exponent <- -(1 / (2 * (1 - rho^2))) * (x1^2 / sig11^2 + x2^2 / sig22^2 - 2 * rho * x1 * x2 / (sig11 * sig22))
    density <- (1 / (2 * pi * sig11 * sig22 * sqrt(1 - rho^2))) * exp(exponent)
    return(density)
  }
  sig11 <- sqrt(sigma[1,1]); sig22 <- sqrt(sigma[2,2]); sig12 <- sigma[1,2]
  rho <- sig12 / (sig11 * sig22)
  # Initialize pis as NA
  pi1 <- NA
  
  # Attempt integration
  tryCatch({
    pi1 <- integrate(function(x1) {
      sapply(x1, function(x1val) {
        integrate(function(x2) bivarnormpdf(x1val, x2, sig11, sig22, rho),
                  lower=0, upper=Inf)$value
      })
    }, lower=0, upper=Inf)$value})
  
  return(c(pi1, 0.25, 0.25, 0.5-pi1))
}


# ## Define a function to compute critical value?
# confint.vpcObj <- function(vpcObj, level=0.95) {
#   z_alph_2 <- -stats::qnorm((1-level)/2)
#   lcl <- vpcObj$vpc - z_alph_2 * vpcObj$stderr
#   ucl <- vpcObj$vpc + z_alph_2 * vpcObj$stderr
#   result <- c(lcl, vpcObj$vpc, ucl)
#   names(result) <- c("lower", "point estimate", "upper")
#   return(result)
# }

####
#formula <- y ~ x + (1 + x|group)
#glmm(formula, nbdata[1,], "nbinom2") 
