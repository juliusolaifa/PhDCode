library(glmmTMB)
library(tweedie)
library(pbapply)

source("GLMMdata.R")
source("GLMMvpc.R")
source("GLMMinference.R")
source("util.R")

glmmsimulatepow <- function(iter, num_per_grp, num_of_grps, beta, sigma, theta, family, 
                            null, alt, power=NULL, vpc_group =1, alpha=0.05) {
  group.sizes <- rep(num_per_grp, num_of_grps)
  x <- rgen01(group.sizes)
  
  if(family == "nbinom2") {
    data <- glmmdata(iter, x, beta, sigma, group.sizes, family, theta=theta)
  } else if(family == "tweedie") {
    data <- glmmdata(iter, x, beta, sigma, group.sizes, family, phi=theta, power=power)
  }
  pvals <- pbsapply(1:nrow(data), function(i) {
    datai <- data[i,]
    lrt(datai, null, alt, family)
  })
  
  powvec <- (pvals < alpha)
  nas <- length(powvec[is.na(powvec)])
  n <- num_per_grp*num_of_grps
  size <- sum(sigma)
  pow <- mean(powvec, na.rm=T)
  result <- c(num_per_grp, num_of_grps, n, size , pow, nas)
  names(result) <- c("n/g", "g", "n", "size", "power", "nas")
  return(result)
}

ss <- expand.grid(5:9,6:10)
ss <- cbind(ss, ss$Var1*ss$Var2)
names(ss) <- c("n/g", "g", "n")
ss <- ss[order(ss$n), ]
set.seed(123)
# Initialize a list to store the results
results_list <- list()

# Loop over each row of ss
for (i in 1:nrow(ss)) {
  # Extract the values from the current row
  n_g <- ss[i, "n/g"]
  g <- ss[i, "g"]
  n <- ss[i, "n"]
  # Call glmmsimulatepow with the current row values
  result <- glmmsimulatepow(iter = 100, 
                            num_per_grp = n_g, 
                            num_of_grps = g, 
                            beta = c(7, 10), 
                            sigma = matrix(c(2, 1, 1, 2), 2), 
                            theta = 2, 
                            family = "nbinom2", 
                            null = y ~ + (1 | group), 
                            alt = y ~ x + (1 + x | group),
                            power = NULL, 
                            vpc_group = 1) 
  
  # Store the result in the results_list
  results_list[[i]] <- result
}
resultdf <- as.data.frame(do.call(rbind, results_list))
