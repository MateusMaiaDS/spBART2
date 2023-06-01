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
n_ <- 200 # Number of observations
# Simulation 1
fried_sim <- mlbench::mlbench.friedman1(n = n_,sd = 0.01)
x <- fried_sim$x[,1:5,drop = FALSE]
x_new <- x
y <- fried_sim$y

# Transforming into data.frame
x <- as.data.frame(x)
x_test <- as.data.frame(x_new)


x_train <- x
y_train <- y
x_test <- x
y_test <- y

continuous_vars <- colnames(x_train)

B_train_arr <- array(data = NA,
                dim = c(nrow(x_train),
                        nrow(knots)+1, # +1 here because is a natural spline
                        ncol(x_train[,continuous_vars, drop = FALSE])))

B_test_arr <- array(data = NA,
               dim = c(nrow(x_test),
                       nrow(knots)+1,  # +1 here because is a natural spline
                       ncol(x_test[,continuous_vars, drop = FALSE])))


# Getting min and max
min_x <- x_min <- apply(as.matrix(x_train),2,min)
max_x <- x_max <- apply(as.matrix(x_train),2,max)


# Getting the internal knots
nIknots <- 10
dif_order <- 2
knots <- apply(x_train,
          2,
          function(x){quantile(x,seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]})



B_train_arr <- array(data = NA,
                dim = c(nrow(x_train),
                        nrow(knots)+1, # +1 here because is a natural spline
                        ncol(x_train[,continuous_vars, drop = FALSE])))

B_test_arr <- array(data = NA,
               dim = c(nrow(x_test),
                       nrow(knots)+1,  # +1 here because is a natural spline
                       ncol(x_test[,continuous_vars, drop = FALSE])))

Z_train_arr <- array(data = NA,
                dim = c(nrow(x_train),
                        nrow(knots)+1-dif_order,
                        ncol(x_train[,continuous_vars, drop = FALSE]))) # correcting the new dimension by P

Z_test_arr <- array(data = NA,
               dim = c(nrow(x_test),
                       nrow(knots)+1-dif_order,
                       ncol(x_test[,continuous_vars, drop = FALSE])))# correcting the new dimension by P


# Creating the natural B-spline for each predictor
for(i in 1:length(continuous_vars)){
B_train_obj <- splines::ns(x = x_train[,continuous_vars[i]],knots = knots[,continuous_vars[i]],
                           intercept = FALSE,
                           Boundary.knots = c(min_x[i],max_x[i]))
B_train_arr[,,i] <- as.matrix(B_train_obj)
B_test_arr[,,i] <- as.matrix(predict(B_train_obj,newx = x_test[,continuous_vars[i]]))
}

# === Directly getting the Pnealised version over the basis function
#see (Eilers, 2010) and look for reference 26 in the text
#=====
if(dif_order!=0){

Z_train_arr <- array(data = NA,
                     dim = c(nrow(x_train),
                             nrow(knots)+1-dif_order,
                             ncol(x_train[,continuous_vars, drop = FALSE]))) # correcting the new dimension by P

Z_test_arr <- array(data = NA,
                    dim = c(nrow(x_test),
                            nrow(knots)+1-dif_order,
                            ncol(x_test[,continuous_vars, drop = FALSE])))# correcting the new dimension by P

D <- D_gen(p = ncol(B_train_arr[,,1]),n_dif = dif_order)

for(i in 1:length(continuous_vars)){
     # IN CASE WE WANT TO USE THE DIFFERENCE PENALISATION DIRECTLY OVER THE
     #BASIS FUNCTION
     Z_train_arr[,,i] <- B_train_arr[,,i]%*%crossprod(D,solve(tcrossprod(D)))
     Z_test_arr[,,i] <- B_test_arr[,,i]%*%crossprod(D,solve(tcrossprod(D)))
}
} else {

# Using original values for B
Z_train_arr <- B_train_arr
Z_test_arr <- B_test_arr
}

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

     for(k in 1:D){
           mean_y[t,k] = mean_y[t,k] + inprod(B[t,1:N_knots,k],beta[1:N_knots,k])

     }

     y[t] ~ dnorm(sum(mean_y[t,1:D]) + gamma, tau)
   }



    # RW prior on beta
    gamma ~ dnorm(0,tau_b_0)

    for(i in 1:N_knots){
         for(j in 1:D){
               beta[i,j] ~ dnorm(0,tau_b[j])
         }
    }

    # Priors on beta values
    tau ~ dgamma(a_tau, d_tau)

    for(k in 1:D){
          tau_b[k] ~ dgamma(0.5 * nu, 0.5 * delta[k] * nu)
          delta[k] ~ dgamma(a_delta, d_delta)
    }

  }
  "

# for(i in 1:length(nu_range)){

# Set up the data
model_data <- list(N = nrow(B_train_arr[,,1]),
              y = as.vector(y_train),
              B = B_train_arr,
              D = ncol(x_train),
              N_knots = ncol(B_train_arr[,,1]),
              # diag_ = diag(nrow = ncol(B_train_arr[,,1])),
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
