rm(list=ls())
library(tidyverse)
library(dbarts)
library(SoftBart)

Rcpp::sourceCpp("src/mpsBART.cpp")
source("R/other_functions.R")
source("R/wrap_bart.R")
source("R/bayesian_simulation.R")
source("R/cv.R")

# Creating the kfold object
n_ <- 100
set.seed(42)

# Simulation 1
sd_ <- 1
n_fold <- 10
# fried_sim <- mlbench::mlbench.friedman1(n = n_,sd = sd_) %>% as.data.frame() %>% .[,c(1:5,11)]
friedman_no_interaction <- function (n, sd = 1)
{
        x <- matrix(runif(1 * n), ncol = 1)
        # y <- 10 * sin(pi * x[, 1] )
        # y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
        y <- 20 * (x[, 1] - 0.5)^2

        if (sd > 0) {
                y <- y + rnorm(n, sd = sd)
        }
        list(x = x, y = y)
}

fried_sim <- friedman_no_interaction(n = n_,sd = sd_) %>% as.data.frame()
k_fold_obj <- k_fold(data = fried_sim,dependent_variable = "y",
                     k_partitions = n_fold,seed = 42,as_data_frame = TRUE)

# Create the data.frame for the RMSE over the test set
model_metrics <- data.frame(n_rep = NA,
                            model = NA,
                            rmse_train = NA,
                            rmse = NA)

for(i in 1:n_fold){

     # Getting the partitions
     x_train <- k_fold_obj[[i]]$x_train
     y_train <- k_fold_obj[[i]]$y_train
     x_test <- k_fold_obj[[i]]$x_test
     y_test <- k_fold_obj[[i]]$y_test


     set.seed(42)
     # Testing the mpsBART
     mpsbart <- rbart(x_train = x_train,y = unlist(c(y_train)),x_test = x_test,
                        n_tree = 1,n_mcmc = 3000,alpha = 0.95,
                        dif_order = 2,
                        beta = 2,nIknots = 30,delta = 1,
                        a_delta = 0.0001,d_delta = 0.0001,nu = 2,a = 0.1,
                        df = 3,sigquant = 0.9,intercept_model = FALSE,stump = TRUE,
                        n_burn = 1000,scale_bool = TRUE)


     # Adding lines in it
     model_metrics <- rbind(model_metrics,data.frame(n_rep = i,
                                                     model = "mpsBART",
                                                     rmse_train = rmse(unlist(y_train),rowMeans(mpsbart$y_hat)),
                                                     rmse = rmse(unlist(c(y_test)),rowMeans(mpsbart$y_hat_test))))

     # Doing the same for BART and
     bart <- dbarts::bart(x.train = x_train,y.train = unlist(y_train),x.test = x_test)

     model_metrics <- rbind(model_metrics,data.frame(n_rep = i,
                                                     model = "BART",
                                                     rmse_train = rmse(unlist(y_train),bart$yhat.train.mean),
                                                     rmse = rmse(unlist(y_test),bart$yhat.test.mean)))

     # Comparing with softBART
     # softBART <- softbart(X = x_train,Y = unlist(y_train),X_test = x_test)
     x_train_soft <- cbind(x_train,x_train)
     x_test_soft <- cbind(x_test,x_test)
     colnames(x_train_soft) <- colnames(x_test_soft) <- c("x1","x2")

     # Returning back to the original with 2 covariates
     # x_train_soft <- x_train
     # x_test_soft <- x_test
     softBART <- softbart(X = x_train_soft,Y = unlist(y_train),X_test = x_test_soft)

     model_metrics <- rbind(model_metrics,data.frame(n_rep = i,
                                                     model = "softBART",
                                                     rmse_train = rmse(unlist(y_train), softBART$y_hat_train_mean),
                                                     rmse = rmse(unlist(y_test),softBART$y_hat_test_mean)))

     # Fitting a regualr spline model
     # mgcv_spline <- gam(y ~ s(x.1,bs = "ps",m = c(2,2))+s(x.2,bs = "ps",m = c(2,2))+s(x.3,bs = "ps",m = c(2,2))+
     #                            s(x.4,bs = "ps",m = c(2,2)) +s(x.5,bs = "ps",m = c(2,2)) +
     #                            te(x.1,x.2), data = cbind(x_train,y_train))
     mgcv_spline <- gam(y ~ s(x,bs = "ps",m = c(2,2))  , data = cbind(x_train,y_train))


     mgcv_pred <- predict(mgcv_spline, newdata = x_test)

     model_metrics <- rbind(model_metrics,data.frame(n_rep = i,
                                                     model = "mgcv",
                                                     rmse_train = rmse(unlist(y_train),mgcv_spline$fitted.values),
                                                     rmse = rmse(unlist(y_test),mgcv_pred)))

     # Fitting a regualr spline model
     # mgcv_spline_no_int <- gam(y ~ s(x.1,bs = "ps",m = c(2,2))+s(x.2,bs = "ps",m = c(2,2))+s(x.3,bs = "ps",m = c(2,2))+
     #                            s(x.4,bs = "ps",m = c(2,2)) +s(x.5,bs = "ps",m = c(2,2)), data = cbind(x_train,y_train))
     # mgcv_spline_no_int <- gam(y ~ s(x,bs = "ps",m = c(2,2)) , data = cbind(x_train,y_train))
     #
     #
     # mgcv_pred_no_int <- predict(mgcv_spline_no_int, newdata = x_test)
     # model_metrics <- rbind(model_metrics,data.frame(n_rep = i,
     #                                                 model = "mgcv_no_intercation",
     #                                                 rmse_train = rmse(unlist(y_train),mgcv_spline_no_int$fitted.values),
     #
     #                                                 rmse = rmse(unlist(y_test),mgcv_pred_no_int)))



}


model_metrics <- model_metrics[complete.cases(model_metrics),]
model_metrics %>%
        ggplot()+
        geom_boxplot(mapping = aes(x = model, y = rmse_train))+
        ggtitle("RMSE over the training set")

model_metrics %>%
ggplot()+
        geom_boxplot(mapping = aes(x = model, y = rmse))+
        ggtitle("RMSE over the test set")


model_metrics %>%
        group_by(model) %>%
        summarise(mean_rmse_test = mean(rmse))

