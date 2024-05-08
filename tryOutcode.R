source("CISim.R")

iter <- 10
num_per_grp <- 5
num_of_grps <- 10
beta <- c(5,7)
sigma <- matrix(c(2,1,1,2), nrow = 2)
theta <- 2
family <- "nbinom2"
fitmodel <- "nbinom2"
formula <- y~x+(1+x|group)
power <- 1.6
vpc_group <- 1

result <- glmmsimulateci(iter = iter, num_per_grp = num_per_grp, 
                         num_of_grps = num_per_grp, beta = beta, 
                         sigma = sigma, theta = theta, family = family, 
                         fitmodel = fitmodel, formula = formula, 
                         power = power, vpc_group = vpc_group)

plot1 <- plotci(result$ci_clean, result$truevpc, result$family, result$fitmodel, type="boundary") 
plot2 <- plotci(result$ciml_clean, result$truevpc, result$family, result$fitmodel, type="non_boundary")

grid.arrange(plot1, plot2, ncol = 2)
