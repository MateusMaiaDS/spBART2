# Header ------------------------------------------------------------------

# P-spline model in JAGS with robust specification of the roughness of the penalty.
# For more details check: "Jullion, Astrid, and Philippe Lambert. Robust
# specification of the roughness penalty prior distribution in spatially
# adaptive Bayesian P-splines models." Computational Statistics & Data Analysis
# 51.5 (2007): 2542-2558."


# Building a JAGS model for splines with the penalized
rm(list=ls())
source("R/other_functions.R")
source("R/wrap_bart.R")
source("R/bayesian_simulation.R")
source("R/cv.R")
library(R2jags)
library(MASS) # Useful for mvrnorm function
library(splines) # Useful for creating the B-spline basis functions
# Maths -------------------------------------------------------------------

# Notation:
# y(t): Response variable at time t, defined on continuous time
# y: vector of all observations
# B: design matrix of spline basis functions
# beta; spline weights

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
set.seed(42)
N <- 200 # Number of observations
x <- sort(runif(N, 0, 10)) # Create some covariate values
nIknots_sim <- 100
knots <- quantile(x,seq(0,1,length.out = nIknots_sim+2))[-c(1,nIknots_sim+2)]
B <- cbind(1,splines::ns(x = x,knots = knots,intercept = FALSE))
tau_b <- 1 # Parameters as above
tau <- 10
beta <- cumsum(c(1, rnorm(ncol(B) - 1, 0, sqrt(tau_b^-1))))
y <- rnorm(N, mean = B %*% beta, sd = sqrt(tau^-1))
# Changing the number of
grid_nIknots <- c(1:10,seq(from = 15, to = 100, by = 5))
iter <- 0
# Iterating over multiple folds-on the test set
n_fold <- 10
n_rep <- 10
rmse_iknot <- matrix(NA, nrow = n_fold*n_rep,ncol = length(grid_nIknots))

for(rep in 1:n_rep){
    data_ <- cbind(x,y)
    k_fold_obj <- k_fold(data = data_,dependent_variable = "y",
                         k_partitions = n_fold,seed = 42+rep)

      for(k in 1:n_fold){
      iter <- iter + 1
      x_train <- k_fold_obj[[k]]$x_train
      y_train <- k_fold_obj[[k]]$y_train
      x_test <- k_fold_obj[[k]]$x_test
      y_test <- k_fold_obj[[k]]$y_test


            for(i in 1:length(grid_nIknots)){
            nIknots <- grid_nIknots[i]
            knots <- quantile(x,seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]
            basis <- splines::ns(x = x_train,knots = knots,intercept = FALSE)
            B_train <- cbind(1,basis)
            B_test <- cbind(1,predict(basis,x_test))
            # plot(x, y)

            # Scaling y
            min_y <- min(y_train)
            max_y <- max(y_train)

            y_scale <- normalize_bart(y = y_train,a = min_y,b = max_y)

            # Setting other parameters
            nsigma <- naive_sigma(x = x_train,y = y_scale)
            df <- 3
            # Calculating tau hyperparam
            a_tau <- df/2
            sigquant <- 0.9

            tau_b_0 <- 4*(2^2)*1
            # Calculating lambda
            qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
            lambda <- (nsigma*nsigma*qchi)/df
            d_tau <- (lambda*df)/2


            # Jags code ---------------------------------------------------------------

            # Jags code to fit the model to the simulated data
            model_code <- "
            model
            {
              # Likelihood
              for (t in 1:N) {
                y[t] ~ dnorm(inprod(B[t,], beta), tau)
              }


              # RW prior on beta
              beta[1] ~ dnorm(0, tau_b_0)
              for (i in 2:N_knots) {
                beta[i] ~ dnorm(beta[i-1], tau_b)
              }

              # Priors on beta values
              tau ~ dgamma(a_tau, d_tau)
              tau_b ~ dgamma(0.5*nu,0.5*delta*nu)
              delta ~ dgamma(a_delta,d_delta)

            }
            "

            # Setting different values for nu
            nu_range <- c(2,5,10,20)
            rough_range <- c(0.00001,0.0001,0.001,0.01)
            # plot(x, y)

            # for(i in 1:length(nu_range)){

              # Set up the data
              model_data <- list(N = nrow(B_train),
                                 y = as.vector(y_train),
                                 B = B_train,
                                 N_knots = ncol(B_train),
                                 a_tau = a_tau,
                                 d_tau = d_tau,
                                 a_delta = 0.0001,
                                 d_delta = 0.0001,
                                 tau_b_0 = tau_b_0,
                                 nu  = 2)

              # Choose the parameters to watch
              model_parameters <- c("beta", "tau", "tau_b","delta")

              # Run the model - can be slow
              model_run <- jags(
                   data = model_data,
                   parameters.to.save = model_parameters,
                   model.file = textConnection(model_code)
              )

              # Simulated results -------------------------------------------------------

              # Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
              # print(model_run)

              # Get the posterior betas and 50% CI
              beta_post <- model_run$BUGSoutput$sims.list$beta
              beta_quantile <- apply(beta_post, 2, quantile, prob = c(0.25, 0.5, 0.75))

              # New prediction
              # B_test <- predict(B_train,x_test)
              y_test_hat <- B_test %*% beta_quantile[2,]
              rmse_iknot[iter, i] <- rmse(x = y_test_hat,y = y_test)
              # Plot the output with uncertainty bands
              # lines(x, B %*% beta, col = alpha("red",0.1)) # True line
              # lines(x, B %*% beta_quantile[2, ], col = i, lty = i) # Predicted line
              # lines(x, B %*% beta_quantile[1, ], col = "blue", lty = 2) # Predicted low
              # lines(x, B %*% beta_quantile[3, ], col = "blue", lty = 2) # Predicted high
              # legend("topleft", c(
              #      "True line",
              #      "Posterior lines (with 50% CI)",
              #      "Data"
              # ),
              # lty = c(1, 1, -1),
              # pch = c(-1, -1, 1),
              # col = c("red", "blue", "black")
              # )
            }
            # Create some new predictions on a grid of new values
            # Needs to be in the same range as the previous values (if not you need to go back to the creation of B above)
            # x_new <- seq(min(x), max(x), length = 1000)
            # nots <- quantile(x,seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]
            # B_new <- cbind(1,splines::ns(x = x_new,knots = knots,intercept = FALSE))
            # plot(x, y)
            # lines(x_new, B_new %*% beta_quantile[2, ], col = "blue") # Beautifully smooth

            }
  }

# Other tasks -------------------------------------------------------------

# Perhaps exercises, or other general remarks

rmse_iknot_summary <- rmse_iknot[rmse_iknot%>% complete.cases(),]
plot(1:28,rmse_iknot_summary %>% colMeans(),xlab = grid_nIknots)
# saveRDS(object = rmse_iknot_summary,file = "R/results/iknot_prior_summary.Rds")
boxplot(rmse_iknot_summary,names = grid_nIknots)
