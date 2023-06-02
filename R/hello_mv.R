rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/spBART.cpp")
source("R/other_functions.R")
source("R/wrap_bart.R")
source("R/bayesian_simulation.R")
n_ <- 500
set.seed(42)

# Simulation 1
sd_ <- 1
fried_sim <- mlbench::mlbench.friedman1(n = n_,sd = sd_)
friedman_no_interaction <- function (n, sd = 1)
{
        # Univariate
        # x <- matrix(runif(1 * n), ncol = 1)
        # y <- 20 * (x[, 1] - 0.5)^2

        x <- matrix(runif(5 * n), ncol = 5)
        y <- 10 * sin(pi * x[, 1]*x[, 2] )
        y <- y + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4] + 5 * x[, 5]

        # x <- matrix(runif(4 * n), ncol = 4)
        # y <- 10 * sin(pi * x[, 1] )
        # y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]

        if (sd > 0) {
                y <- y + rnorm(n, sd = sd)
        }
        list(x = x, y = y)
}

fried_sim <- friedman_no_interaction(n = n_,sd = sd_)
fried_sim_new_sample <- friedman_no_interaction(n = n_,sd = sd_)

x <- fried_sim$x[,,drop = FALSE]
x_new <- fried_sim_new_sample$x
y <- fried_sim$y

# Transforming into data.frame
x <- as.data.frame(x)
x_test <- as.data.frame(x_new)


# Testing the mpsBART
bart_test <- rbart(x_train = x,y = unlist(c(y)),x_test = x_test,
                   n_tree = 50,n_mcmc = 2500,alpha = 0.95,
                   dif_order = 0,
                   beta = 2,nIknots = 5,delta = 1,
                   a_delta = 0.0001,d_delta = 0.0001,nu = 2,
                   df = 3,sigquant = 0.9,a = 0,
                   n_burn = 500,scale_bool = TRUE)


# Running BART
bartmod <- dbarts::bart(x.train = x,y.train = unlist(c(y)),ntree = 200,x.test = x_test,keeptrees = TRUE)

# Convergence plots
par(mfrow = c(1,2))
plot(bart_test$tau_post,type = "l", main = expression(tau),ylab=  "")
plot(bartmod$sigma^-2, type = "l", main = paste0("BART: ",expression(tau)),ylab=  "")

par(mfrow = c(1,2))
plot(bart_test$y_hat %>% rowMeans(),y, main = 'mpsBART', xlab = "mpsBART pred", ylab = "y")
plot(bartmod$yhat.train.mean,y, main = "BART", xlab = "BART pred", ylab = "y")

# Comparing on the test set
pred_bart <- colMeans(predict(bartmod,fried_sim_new_sample$x))

# Storing the results
rmse(x = fried_sim_new_sample$y,y = bartmod$yhat.test.mean)
rmse(x = fried_sim_new_sample$y,y = rowMeans(bart_test$y_hat_test))


# Storing the results
rmse(x = fried_sim$y,y = bartmod$yhat.train.mean)
rmse(x = fried_sim$y,y = rowMeans(bart_test$y_hat))

par(mfrow=c(1,1))
plot(bartmod$yhat.test.mean,rowMeans(bart_test$y_hat_test))
plot(bartmod$yhat.train.mean,rowMeans(bart_test$y_hat))

# col_zero_two <- apply(x_train[,-1],1,sum)
# plot(x$V1,y)
# points(x_test$V1,bartmod$yhat.test.mean,col = "blue")
# points(x_test$V1,bart_test$y_hat_test %>% rowMeans(),col = "orange")
