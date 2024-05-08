library(glmmTMB)
library(tweedie)
library(pbapply)

source("GLMMdata.R")
source("GLMMvpc.R")
source("GLMMinference.R")
source("util.R")

glmmsimulatealp <- function(iter, num_per_grp, num_of_grps, beta, sigma, theta, family, 
                            formula, alt, power=NULL, vpc_group =1, alpha=0.05) {
  group.sizes <- rep(num_per_grp, num_of_grps)
  x <- rgen01(group.sizes)
  
  if(family == "nbinom2") {
    data <- glmmdata(iter, x, beta, sigma, group.sizes, family, theta=theta)
  } else if(family == "tweedie") {
    data <- glmmdata(iter, x, beta, sigma, group.sizes, family, phi=theta, power=power)
  }
  alp <- pbsapply(1:nrow(data), function(i) {
    print(i)
    datai <- data[i,]
    fit <- glmm(formula, datai, family, computepis = F)
    lrt(fit, alt)
  })
  noofNA <- length(alp[is.na(alp)])
  result <- c(mean(alp<alpha, na.rm=T), nrow(data)-noofNA, nrow(data))
  names(result) <- c("alpha", "no of features", "no of Iterations")
  return(result)
}
set.seed(123)
out <- glmmsimulatealp(iter=1000, num_per_grp=10, num_of_grps=20, beta=c(7,10), 
                       sigma=matrix(c(2,0,0,0),2), theta=0.5, family="nbinom2", 
                       formula=y ~ x + (1 |group), alt = y ~ x + (1 + x|group),
                       power=NULL, vpc_group =1) 
