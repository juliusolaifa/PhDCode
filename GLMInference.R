####### CONFIDENCE INTERVAL  #######
library(tmvtnorm)
#Projection Matrix
computeProjectionMatrix <- function(basisMatrix, Itheta) {
  I <- diag(nrow(Itheta))
  I - Itheta %*% basisMatrix %*% MASS::ginv(t(basisMatrix) %*% Itheta %*% basisMatrix) %*% t(basisMatrix)
}


monte_carlo_quantile <- function(pis, grad, Itheta, alpha, N) {
  
  #Basis Matrices as defined by Self & Liang 1987
    n <- nrow(Itheta)
    if(length(pis)==4) {
        basis <- list(
          matrix(0, nrow = n, ncol = n), # B1
          diag(c(rep(0, 3),1,rep(0,n-4))),     # B2
          diag(c(rep(0, 5),1,rep(0,n-6))),  # B3
          diag(c(rep(0, 3),1,0,1,rep(0,n-6)))   # B4
        )
    } else if(length(pis) == 2) {
      basis <- list(
          matrix(0, nrow = n, ncol = n), # B1
          diag(c(rep(0, 3),1,rep(0,n-4)))     # B2
      )
    }
    
    #Generate a large sample
    r <- runif(N)
    
    #This trick saves from using multiple if statements
    pis <- c(0,cumsum(pis))
    indices <- findInterval(r, pis)
  
    W <- vector()
   
    # Based on the realization of r, determine which of the mixture
    # you will simulate from
    for(j in seq_along(indices)) {
      index <- indices[j]
      B <- basis[[index]]
      P <- computeProjectionMatrix(B,Itheta)
      # no of parameters bounded
      q <- sum(diag(B) == 1)
      mean_vec <- rep(0, n - q)
      sigma_mat <- (P %*% Itheta %*% t(P))
      
  
      # which parameter is bounded, remove from the information matrix
      cond <- which(diag(B) == 1)
      if(length(cond)>0)
        sigma_mat <- sigma_mat[-cond, -cond]
      
      #(n-q) variate
      lower_bounds <- rep(-Inf, n - q)
      upper_bounds <- rep(Inf, n - q)
      w <- as.vector(rtmvnorm(n = 1, mean = mean_vec, sigma = sigma_mat,
                    lower = lower_bounds, upper = upper_bounds))
      if(q == 1){
        w <- append(w, 0, after = cond-1)
      }else if(q == 2){
        w <- append(w, 0, after = 3)
        w <- append(w, 0, after = 5)
      }
      w <- t(grad) %*% w
      W <- c(W, w)
    }
    quantile(W, probs=c(alpha/2, 1-alpha/2))
}

c.i_vpc <- function(vpcObj, type, alpha = 0.05, N = 1000) {
  if(is.null(vpcObj)) {
    c.i <- c(NA, NA, NA)
  }else{
    if(type=="non_boundary") {
      sterr <- sqrt((vpcObj$gradient %*% vpcObj$Itheta %*% vpcObj$gradient)/vpcObj$n)
      crit_val <- qnorm(1-alpha/2)
      c.i <- c(vpcObj$vpc - crit_val*sterr, vpcObj$vpc, vpcObj$vpc + crit_val*sterr)
    } else if(type=="boundary"){
      w <- monte_carlo_quantile(vpcObj$pis, vpcObj$gradient, vpcObj$Itheta,
                                alpha, N)/sqrt(vpcObj$n)
      c.i <- c(vpcObj$vpc - w[2], vpcObj$vpc, vpcObj$vpc-w[1])
    }
  }
  names(c.i) <- c("lcl", "est", "ucl")
  c.i
}

####### HYPOTHESIS TESTING #######
lrt <- function(data, H0=NULL, H1=NULL, family) {
  glmmH1 <- glmmTMB(H1, data, family)
  glmmH0 <- glmmTMB(H0, data, family)  
  if(is.na(logLik(glmmH1)) && is.na(logLik(glmmH0))) {
    return(NA)
  }else if(is.na(logLik(glmmH1)) && !is.na(logLik(glmmH0))){
    return(1)
  }

  test.stat <- 2*(as.numeric(logLik(glmmH1))- as.numeric(logLik(glmmH0)))
   pval <- as.numeric(0.5*pchisq(test.stat, df = 1, lower.tail = FALSE) +
                           0.5*pchisq(test.stat, df = 2, lower.tail = FALSE))
  pval
}
