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
n_ <- 250
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
                            rmse = NA,
                            n_tree = NA)
number_trees <- 1:10

# Getting y pred
# y_pred_train_tree <- matrix(NA, nrow = n_, ncol = length(number_trees))
y_pred_test_tree <- matrix(NA, nrow = n_, ncol = length(number_trees))

for(n_t in 1:length(number_trees)){
     y_pred_train <- numeric()
     y_pred_test <- numeric()

     for(i in 1:n_fold){

          # Getting the partitions
          x_train <- k_fold_obj[[i]]$x_train
          y_train <- k_fold_obj[[i]]$y_train
          x_test <- k_fold_obj[[i]]$x_test
          y_test <- k_fold_obj[[i]]$y_test


          set.seed(42)
          # Testing the mpsBART
          mpsbart <- rbart(x_train = x_train,y = unlist(c(y_train)),x_test = x_test,
                           n_tree = number_trees[n_t],n_mcmc = 2500,alpha = 0.95,
                           dif_order = 2,
                           beta = 2,nIknots = 30,delta = 1,
                           a_delta = 0.0001,d_delta = 0.0001,nu = 2,
                           df = 3,sigquant = 0.9,intercept_model = FALSE,a = number_trees[n_t],
                           n_burn = 500,scale_bool = TRUE)


          # Adding lines in it
          model_metrics <- rbind(model_metrics,data.frame(n_rep = i,
                                                          model = "mpsBART",
                                                          n_tree = number_trees[n_t],
                                                          rmse_train = rmse(unlist(y_train),rowMeans(mpsbart$y_hat)),
                                                          rmse = rmse(unlist(c(y_test)),rowMeans(mpsbart$y_hat_test))))

          y_pred_train <- c(y_pred_train,rowMeans(mpsbart$y_hat))
          y_pred_test <- c(y_pred_test,rowMeans(mpsbart$y_hat_test))


     }
     # y_pred_train_tree[,n_t] <- y_pred_train
     y_pred_test_tree[,n_t] <- y_pred_test

     print(paste(" Tree number iteration value: ", number_trees[n_t]))
}


model_metrics <- model_metrics[complete.cases(model_metrics),] %>%
                    mutate(n_tree = as.factor(n_tree))
model_metrics %>%
     ggplot()+
     geom_boxplot(mapping = aes(x = n_tree, y = rmse))+
     ggtitle("RMSE over the test set (mpsBART)")

model_metrics %>%
     ggplot()+
     geom_boxplot(mapping = aes(x = n_tree, y = rmse_train))+
     ggtitle("RMSE over the training set (mpsBART)")

model_metrics %>%
     group_by(model) %>%
     summarise(mean_rmse_test = mean(rmse))
par(mfrow=c(1,1))
plot(unlist(purrr::map(k_fold_obj, ~.x$x_test)),unlist(purrr::map(k_fold_obj, ~.x$y_test)), xlab = "x", ylab = "y")
for(i in 1:length(number_trees)){
     if(number_trees[i]<=5){
          col <- "black"
     } else {
          col <- "red"
     }
     points(unlist(purrr::map(k_fold_obj, ~.x$x_test)),y_pred_test_tree[,i], pch = 20, col = alpha(i,0.5))
}
