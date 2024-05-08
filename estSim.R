library(pbapply)

source("GLMMdata.R")
source("GLMMvpc.R")
source("GLMMinference.R")
source("util.R")

glmmsimulateest <- function(iter, num_per_grp, num_of_grps, beta, sigma, theta, family, 
                           fitmodel, formula, power=NULL) {
  group.sizes <- rep(num_per_grp, num_of_grps)
  x <- rgen01(group.sizes)
  
  if(family == "nbinom2") {
    data <- glmmdata(iter, x, beta, sigma, group.sizes, family, theta=theta)
  } else if(family == "tweedie") {
    data <- glmmdata(iter, x, beta, sigma, group.sizes, family, phi=theta, power=power)
  }
  truevpc0 <- switch(family,
                    "nbinom2" = vpc_nbinom2(beta[1], beta[2], theta, sigma[1,1], sigma[1,2], 
                                            sigma[2,2], x=0),
                    "tweedie" = vpc_tweedie(beta[1], beta[2], phi=theta, sigma[1,1], sigma[1,2], 
                                            sigma[2,2],power, x=0)
  )
  truevpc1 <- switch(family,
                    "nbinom2" = vpc_nbinom2(beta[1], beta[2], theta, sigma[1,1], sigma[1,2], 
                                            sigma[2,2], x=1),
                    "tweedie" = vpc_tweedie(beta[1], beta[2], phi=theta, sigma[1,1], sigma[1,2], 
                                            sigma[2,2],power, x=1)
  )
  
  fit <- as.data.frame(t(pbsapply(1:nrow(data), function(i) {
    print(i)
    datai <- data[i,]
    fit <- glmm(formula, datai, fitmodel)
    vpc0 <- vpc(fit, x=0)$vpc
    vpc1 <- vpc(fit, x=1)$vpc
    if (is.null(vpc0) || is.null(vpc1)) {
      return(list(vpc0=NA, vpc1=NA))
    } else {
      return(list(vpc0=vpc0, vpc1=vpc1))
    }
  })))
  

  
  fit$vpc0 <- as.numeric(fit$vpc0)
  fit$vpc1 <- as.numeric(fit$vpc1)
  fit$truevpc0 <- as.numeric(truevpc0)
  fit$truevpc1 <- as.numeric(truevpc1)
  fit$bias0 <- fit$vpc0-fit$truevpc0
  fit$bias1 <- fit$vpc1-fit$truevpc1
# 
#   boxplot(fit$vpc0, fit$vpc1,
#           names = c("VPC0", "VPC1"),
#           main = "VPC Boxplots",
#           ylab = "VPC Values",
#           col = c("pink", "lightblue"))
  
  # out <- cbind(fit$vpc0-fit$truevpc0, fit$vpc1-fit$truevpc1)
  # colnames(out) <- c("Bias0", "Bias1")
  return(fit)
}

plotest <- function(biases) {
  boxplot(biases)
}

# Create the DataFrame in R
parameters <- data.frame(
  b0 =    c(5, 5, 2, 2, 2, 2, 2, 2, 2, 2),
  b1 =    c(4, 4, 3, 3, 3, 3, 3, 3, 3, 3),
  sig11 = c(4, 5, 3, 2, 4, 3, 3, 2, 2, 2),
  sig12 = c(2, 1, 2, 1, 1, 1, 2, 1, 1, 1),
  sig22 = c(5, 5, 3, 2, 2, 2, 2, 2, 2, 2),
  theta = c(0.1, 0.25, 0.45, 0.7, 1, 1.5, 2.4, 4, 9, 45)
)
#theta = c(0.10, 0.20, 0.30, 0.40, 0.75, 1.00, 1.50, 2.00, 5.00, 15.00)
#  VPC = c(0.861641, 0.758280, 0.677061, 0.611558, 0.456859, 0.386944, 0.296266, 0.240019, 0.112205, 0.040433),
#VPC_rounded = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0)

# b0 <- 2
# b1 <- 3
# sig11 <- c(2,3)
# sig12 <- c(-1,0,1)
# sig22 <- c(1,2)
# theta <- c(0.1,0.25,0.5,1,2,4,6,20)
# 
# parameters <- expand.grid(b0=b0, b1=b1, sig11=sig11, sig12=sig12,
#                           sig22=sig22,theta=theta)

set.seed(123)
# result <- glmmsimulateest(iter = 10, num_per_grp = 10, num_of_grps = 20, beta = c(5,4),
#                sigma = matrix(c(3,2,2,2),2), theta=45, family="nbinom2",
#                fitmodel="nbinom2", formula=y~x+(1+x|group), power=1.6)


apply_results <- apply(parameters, 1, function(row) {
  # Extract the parameters from the row
  b0 <- row["b0"]
  b1 <- row["b1"]
  sig11 <- row["sig11"]
  sig12 <- row["sig12"]
  sig22 <- row["sig22"]
  theta <- row["theta"]

  # Call glmmsimulateest with the extracted parameters
  result <- glmmsimulateest(iter = 100, num_per_grp = 10, num_of_grps = 20, beta = c(b0, b1),
                            sigma = matrix(c(sig11, sig12, sig12, sig22), 2), theta = theta,
                            family = "nbinom2", fitmodel = "tweedie", formula = y ~ x + (1 + x | group),
                            power = 1.6)
  result <- as.data.frame(result)

  # Return the result
  return(result)
})

# final_result <- do.call(rbind, apply_results)
# 
# plot(final_result$truevpc0, final_result$bias0, xlim=c(0,1), ylim=c(-1,1), xlab="VPC", ylab="Bias")

# Set up the layout for the plots
par(mfrow=c(2,5))

# Iterate through each element in apply_results
for (i in 10:1) {
  # Extract vpc0 and vpc1 from apply_results[[i]]
  vpc0 <- apply_results[[i]]$vpc0
  vpc1 <- apply_results[[i]]$vpc1
  
  # Combine vpc0 and vpc1 into a single data frame
  vpc_data <- data.frame(
    vpc0 = vpc0,
    vpc1 = vpc1
  )
  
  # Create a box plot for vpc0 and vpc1
  boxplot(vpc_data, 
          names = c("vpc0", "vpc1"),
          main = paste("Boxplot of vpc0 and vpc1 for pair", i),
          ylab = "VPC Values",
          col = c("pink", "lightblue"))
  
  # Add true VPC values
  points(x = 1:2, y = c(apply_results[[i]]$truevpc0[1], apply_results[[i]]$truevpc1[1]), col = c("red", "blue"), pch = 19)
}

