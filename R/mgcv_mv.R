rm(list=ls())
library(tidyverse)
library(mgcv)
library(mlbench)
library(dbarts)

source("R/other_functions.R")
n_ <- 100
set.seed(42)

no_interaction_friedman <- function (n, sd = 1)
{
     x <- matrix(runif(10 * n), ncol = 4)
     y <- 10 * sin(pi * x[, 1] )
     y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
     if (sd > 0) {
          y <- y + rnorm(n, sd = sd)
     }
     list(x = x, y = y)
}

# Simulation 1
fried_sim <- no_interaction_friedman(n = n_,sd = 1)
x <- fried_sim$x
x_new <- x
y <- fried_sim$y

# Transforming into data.frame
x <- as.data.frame(x)
x_test <- as.data.frame(x_new)

data_ <- data.frame(x,y)

gam_mod <- gam(y ~ s(V1)+s(V2)+s(V3)+s(V4),data = data_)
pred_gam <- predict(gam_mod,x)

# Getting value of sig2
gam_mod$sig2


# Comparing with BART
bartmod <- bart(x.train = x,y.train = y)

# Comparing predictions
par(mfrow = c(1,2))
plot(pred_gam,y, main = 'GAM', xlab = "GAM pred", ylab = "y")
plot(bartmod$yhat.train.mean,y, main = "BART", xlab = "BART pred", ylab = "y")


