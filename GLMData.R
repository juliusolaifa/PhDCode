source("util.R")
library(glmmTMB)
library(tweedie)


glmmdata <- function(iter.size, beta, sigma.u, group.sizes, 
                     family, x = NULL, theta = NULL, 
                     phi = NULL, power = NULL) {
  
  #Create the fixed effect design matrix X
  if(!is.null(x)) {
      x <- as.matrix(x)
      if(nrow(x) != sum(group.sizes)) {
        stop("number of rows of x must be equal to sum of cluster sizes")
      }
      X <- stats::model.matrix(~x)
  }
  
  if(is.matrix(sigma.u) && is.null(x)) {
    stop("x ca not be null for a slope model")
  }
  
  # Obtain the number of groups
  n_groups <- length(group.sizes)
  
  # Generate a the random effect, intercept or slope model
  if(is.matrix(sigma.u)) {
      n_rand_par <- ncol(sigma.u)
      u <- MASS::mvrnorm(n=n_groups, mu=rep(0,n_rand_par), Sigma = sigma.u)
  } else {
      u <- rnorm(n=n_groups, mean = 0, sd=sigma.u)
  }
  
  # Create a variable group as factor from the group.size vector
  group_names <- paste0("group", 1:length(group.sizes))
  group_assignments <- rep(x=group_names, times=group.sizes)
  group_assignments <- factor(group_assignments)

  
  # Create the random effect design matrix Z
  z0 <- stats::model.matrix(~ 0 + group_assignments)
  if(!is.null(x) && is.matrix(u)) {
        z <- do.call(cbind, lapply(1:ncol(u), function(i) z0 * X[, i])) 
  } else{
        z <- z0
  }

  # Create a matrix to hold the response variable
  total_sample_size <- sum(group.sizes)
  y_matrix <- matrix(nrow=iter.size, ncol=total_sample_size)

  # Obtain the conditional mean 
  if(!is.null(x)){
      mu <- exp(X%*%beta + z%*%c(u))
  }else{
      mu <- exp(beta + z%*%c(u))
  }
  
  # Split the conditional mean into the groups
  cond_means <- split(x =mu, f = rep(1:length(group.sizes), group.sizes))
  
  # Validate that the nuisance parameter(s) of the conditional distribution of y
  # is/are supplied
  if(family == "nbinom2" && is.null(theta)) {
    stop("theta can not be NULL when family is nbinom2")
  }else if(family == "tweedie" && (is.null(phi) || is.null(power)) ){
    stop("phi or power can not be NULL when family is tweedie")
  }

  # Generate data
  generate_y <- switch(family,
     nbinom2={
       function() mapply(MASS::rnegbin, n=group.sizes, 
                                mu=cond_means,theta=theta)
     },
     tweedie={
       function() mapply(tweedie::rtweedie, n=group.sizes, 
                                mu=cond_means, phi=phi, power=power)
     },
     stop("data generation not implemented for family: ", family)
  )
  
  
  # Store the data in the y_matrix, Can we remove the loop?
  for(i in 1:iter.size) {
    y_matrix[i, ] <- c(generate_y())
  }
  result <- structure(list("x" = x, "z" = z, "y" = y_matrix, "family" = family, 
                           "group" = group_assignments), class = "glmmData")
  result
}
# indexing method for `heritData` class
"[.glmmData" <- function(dataObj, i, ...) {
  if(!is.null(dataObj$x)) {
    data <- data.frame(x = dataObj$x, group = dataObj$group, y = dataObj$y[i,])
  } else {
    data <- data.frame(group = dataObj$group, y = dataObj$y[i,])
  }
}

print.glmmData <- function(dataObj) {
  cat("Family:", dataObj$family, "\n")
  
  if(is.null(data$x)) {
      datamat <- rbind(t(dataObj$x), dataObj$y)
      nx <- ncol(dataObj$x)
      rownames(datamat) <- c(paste0('x', 1:nx), paste0('row', 1:nrow(dataObj$y)))
  } else{
      datamat <-  dataObj$y
      rownames(datamat) <- paste0('row', 1:nrow(dataObj$y))
  }
  colnames(datamat) <- dataObj$group
  
  print(datamat)
}

summary.glmmData <- function(dataObj) {
  cat("Summary of glmmData object\n")
  cat("----------------------------\n")
  cat("Family:", dataObj$family, "\n\n")
  cat("Dimensions of data:", dim(dataObj)[1], "rows,", dim(dataObj)[2], "columns\n\n")
  cat("Number of groups:", length(unique(dataObj$group)), "\n")
  cat("Group sizes:\n")
  print(table(dataObj$group))
  invisible(dataObj)
}

dim.glmmData <- function(datObj) {
  dim(datObj$y)
}

nrow.glmmData <- function(dataObj) {
  dim(dataObj)[1]
}

head.glmmData <- function(dataObj, n = 5L, ...) {
  head_y <- utils::head(dataObj$y, n)
  structure(list(x = dataObj$x, y = head_y, family = dataObj$family, 
                 group = dataObj$group), class = "glmmData")
}

tail.glmmData <- function(dataObj, n = 5L, ...) {
  tail_y <- utils::tail(dataObj$y, n)
  structure(list(x = dataObj$x, y = tail_y, family = dataObj$family, 
                 group = dataObj$group), class = "glmmData")
}

# Example
#x <- c(23,24,25,26,29,30,18,22,21,16,18,16,21)
group.sizes <- rep(10,20)
x <- rgen01(group.sizes)
beta <- c(7,9)
sigma.u <- matrix(c(2,1,1,2),2)
theta = 2
phi = 2
power = 1.6
nu = 3
lambda2 = 0.7
iter =1000
nbdata <- glmmdata(iter,beta, sigma.u, group.sizes, "nbinom2", x= x, theta=theta)
#twdata <- glmmdata(iter, beta, sigma.u, group.sizes, "tweedie", x = x, phi=phi, power=power)
# data3 <- glmmdata(iter, x, beta, sigma.u, group.sizes, "compois", nu=nu)
