set.seed(123)
library(pbapply)
library(ggplot2)
library(gridExtra)

source("GLMMdata.R")
source("GLMMvpc.R")
source("GLMMinference.R")
source("util.R")

glmmsimulateci <- function(iter, num_per_grp, num_of_grps, beta, sigma, theta, family, 
                           fitmodel, formula, power=NULL, vpc_group) {
  group.sizes <- rep(num_per_grp, num_of_grps)
  x <- rgen01(group.sizes)
  
  if(family == "nbinom2") {
    data <- glmmdata(iter, beta, sigma, group.sizes, family, x=x, theta=theta)
  } else if(family == "tweedie") {
    data <- glmmdata(iter, beta, sigma, group.sizes, family, x=x, phi=theta, power=power)
  }
  
  if(is.matrix(sigma) && dim(sigma)[1] == 2){
    truevpc <- switch(family,
                      "nbinom2" = vpc_nbinom2(beta[1], beta[2], theta, sigma[1,1], sigma[1,2],
                                              sigma[2,2], x=1),
                      "tweedie" = vpc_tweedie(beta[1], beta[2], phi=theta, sigma[1,1], sigma[1,2],
                                              sigma[2,2],power, x=1)
    )
  }else if(length(sigma) == 1) {
    truevpc <- switch(family,
                      "nbinom2" = vpc_nbinom2(beta[1], beta[2], theta, sigma, 0,
                                              0, x=1),
                      "tweedie" = vpc_tweedie(beta[1], beta[2], phi=theta, sigma, 
                                              0, 0, power, x=1)
    )
  }
  
  ci <- as.data.frame(t(pbsapply(1:nrow(data), function(i) {
    datai <- data[i,]
    fit <- glmm(formula, datai, fitmodel)
    vpcc <- vpc(fit, x=vpc_group)
    c.i_vpc(vpcc, type="boundary")
  })))
  
  print(ci)
  
  ci_clean <- ci[complete.cases(ci) & !(ci$lcl < 0) & !(ci$ucl > 1),]
  
  ciml <- as.data.frame(t(pbsapply(1:nrow(data), function(i) {
    datai <- data[i,]
    fit <- glmm(formula, datai, fitmodel)
    vpcc <- vpc(fit, x=vpc_group)
    c.i_vpc(vpcc, type="non_boundary")
  })))
  print(ciml)
  ciml_clean <- ciml[complete.cases(ciml) & !(ciml$lcl < 0) & !(ciml$ucl > 1),]
  
  return(list("ci_clean"=ci_clean, "ciml_clean"=ciml_clean, "truevpc" = truevpc, 
              "family" = family, "fitmodel" = fitmodel))
}

familyname <- function(name) {
  if(name == "nbinom2") "Negative Binomial"
  else "Tweedie"
}

plotci <- function(ci, truevpc,family,fitmodel, type) {
  cov_prob <- sum(sapply(1:nrow(ci), function(i) {
    ci[i,1] < truevpc && truevpc < ci[i,3]
  })) / nrow(ci)  # Adjusted to divide by number of rows for proper probability
  
  cov_width <- mean(sapply(1:nrow(ci), function(i) {
    ci[i,3] - ci[i,1]
  }))
  
  truevpc <- as.numeric(truevpc)
  rtruevpc <- round(truevpc, 2)
  ci$color_indicator <- ifelse(ci$lcl <= truevpc & ci$ucl >= truevpc, "blue", "red")
  
  p <- ggplot(ci, aes(x = 1:nrow(ci))) +
    geom_segment(aes(xend = 1:nrow(ci), y = lcl, yend = ucl, color = color_indicator), linewidth = 0.5) +
    geom_point(aes(y = lcl), color = "black", size = 1) +  # Point at the lower confidence limit
    geom_point(aes(y = ucl), color = "black", size = 1) +  # Point at the upper confidence limit
    geom_point(aes(y = est), color = "lightgreen", size = 2) +
    geom_hline(yintercept = min(ci$est), color = "grey", linetype = "solid", linewidth = 0.5) +
    geom_hline(yintercept = max(ci$est), color = "grey", linetype = "solid", linewidth = 0.5) +
    geom_hline(yintercept = truevpc, color = "black", linetype = "solid", linewidth = 0.5) +
    scale_color_identity() +
    scale_y_continuous(
      breaks = c(min(ci$lcl), min(ci$est), truevpc, max(ci$est), max(ci$ucl)),  # Custom breaks including min/max estimates
      labels = c(paste("Min LCL:", round(min(ci$lcl), 2)),
                 paste("Min Est:", round(min(ci$est), 2)),
                 paste("True VPC:", round(truevpc, 2)),
                 paste("Max Est:", round(max(ci$est), 2)),
                 paste("Max UCL:", round(max(ci$ucl), 2))
      ),  # Custom labels for the breaks
      limits = c(min(ci$lcl), max(ci$ucl))  # Adjust limits to fit min and max of CI
    ) +
    labs(x = paste("Confidence Interval", type),
         y = "Variance Partition Coefficient",
         title = paste("Data Generating Model:", familyname(family), "\n\t\tFitted Model:", familyname(fitmodel)),
         subtitle = paste("95% Confidence Interval \nCoverage probability =", round(cov_prob,4), "\nAverage Coverage width =", round(cov_width, 4))) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(), legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}
