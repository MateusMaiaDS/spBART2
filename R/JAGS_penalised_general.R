# Header ------------------------------------------------------------------

# Sum of P-spline model in JAGS with robust specification of the roughness of the penalty.
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
library(tidyverse)

# Maths -------------------------------------------------------------------

# Notation:
# y(x): Response variable at variable, defined on continuous x
# y: vector of all observations
# B: design matrix of spline basis functions
# beta; spline weights
# n_sp: number of splines
# P: Penalised difference matrix


# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
set.seed(42)
# Creating the kfold object
n_ <- 100
set.seed(42)

# ============
# Simulation 1
# ============
sd_ <- 1

friedman_univariate<- function (n, sd = 1)
{
     x <- matrix(runif(1 * n), ncol = 1)
     y <- 20 * (x[, 1] - 0.5)^2

     if (sd > 0) {
          y <- y + rnorm(n, sd = sd)
     }
     list(x = x, y = y)
}

fried_sim <- friedman_univariate(n = n_,sd = sd_) %>% as.data.frame() %>% arrange(x)

x <- fried_sim$x
y <- fried_sim$y

# Setting
nIknots <- 20
# knots <- quantile(x,seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]

# To understand the Boundary Knots see rejoinder sectin from Eilers (1996)

tau_b <- 1 # Parameters as aboveau_b^-1))))
min_x <- min(x)
max_x <- max(x)
ndx <- nIknots
bdeg <- 3
dx <- (max_x-min_x)/ndx
knots <- seq(from = min_x-bdeg*dx,to = max_x+bdeg*dx, by = dx)
basis <- splines::spline.des(knots = knots,x = x,ord = bdeg+1,derivs = x*0,outer.ok = FALSE)$design
B_train <- basis


# Setting x and y to be x_train and x_test
x_train <- x_test <- x
y_train <- y

# Scaling y
min_y <- min(y_train)
max_y <- max(y_train)

y_scale <- normalize_bart(y = y_train,a = min_y,b = max_y)

# Setting the sum elements
n_ps <- 10
ones <- matrix(1,nrow = n_ps)

# Setting other parameters
nsigma <- naive_sigma(x = x_train,y = y_scale)
df <- 3
# Calculating tau hyperparam
a_tau <- df/2
sigquant <- 0.9
kappa <- 2
tau_b_0 <- 4*(kappa^2)*n_ps
# Calculating lambda
qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
lambda <- (nsigma*nsigma*qchi)/df
d_tau <- (lambda*df)/2


# Calculating penalisation matrix
n_diffs <- 2
if(n_diffs==0){
     D <- diag(1,nrow = ncol(B_train))
} else {
     D <- diff(diag(ncol(B_train)), diff = n_diffs )
}

P <- crossprod(D)
P_0 <- matrix(0 , nrow = nrow(P), ncol = ncol(P))

if(n_diffs >1 ){
  diag(P_0[1:n_diffs,1:n_diffs]) <- 1
} else{
  P_0[1,1] <- 1
}


# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code <- "
model {
  # Likelihood
  y[1:N] ~ dmnorm((B%*%betas)%*%ones+rep(beta_intercept,N),tau*diag_)

  # RW prior on beta
  beta_intercept ~ dnorm(0,tau_b_0)

  for(i in 1:n_ps){
     betas[1:N_knots,i] ~ dmnorm(rep(0,N_knots),tau_b*(a*P_0 + P))
  }

  # Priors on beta values
  tau_b ~ dgamma(0.5 * nu, 0.5 * delta * nu)
  tau ~ dgamma(a_tau, d_tau)
  delta ~ dgamma(a_delta, d_delta)
}"



# Set up the data
model_data <- list(N = nrow(B_train),
                   n_ps = n_ps,
                   diag_ = diag(1,nrow = nrow(B_train)),
                   y = as.vector(y_train),
                   B = B_train,
                   P = P,
                   P_0 = P_0,
                   ones = ones,
                   N_knots = ncol(B_train),
                   a_tau = a_tau,
                   d_tau = d_tau,
                   a_delta = 0.0001,
                   d_delta = 0.0001,
                   tau_b_0 = tau_b_0,
                   nu  = 2,
                   a = 0.0001)

# Choose the parameters to watch
model_parameters <- c("betas", "tau", "tau_b","delta","beta_intercept")

# Run the model - can be slow
model_run <- jags(
     data = model_data,
     parameters.to.save = model_parameters,
     model.file = textConnection(model_code)
)

# Simulated results -------------------------------------------------------
saveRDS(model_run,file = paste0("R/JAGS_exp/model_run_",n_diffs,"_dif.Rds"))
# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
# print(model_run)

# Get the posterior betas and 50% CI
beta_post <- model_run$BUGSoutput$sims.list$betas
beta_mean <- apply(beta_post, c(2,3), mean)
beta_low <- apply(beta_post, c(2,3), function(x)quantile(x,prob = c(0.025)))
beta_up <- apply(beta_post, c(2,3), function(x)quantile(x,prob = c(0.975)))
each_tree_pred <- matrix(0,nrow  = nrow(B_train), ncol = n_ps)

# New prediction
y_train_hat <- (B_train %*% beta_mean) %*% ones + matrix(model_run$BUGSoutput$mean$beta_intercept, nrow = nrow(B_train))
y_train_hat_low <- (B_train %*% beta_low) %*% ones + matrix(model_run$BUGSoutput$mean$beta_intercept, nrow = nrow(B_train))
y_train_hat_up <- (B_train %*% beta_up) %*% ones + matrix(model_run$BUGSoutput$mean$beta_intercept, nrow = nrow(B_train))

# rmse_iknot[iter, i] <- rmse(x = y_test_hat,y = y_test)
# Plot the output with uncertainty bands
par(mfrow=c(1,1))
plot(x,y, main = c("Sum of splines"))

for(i in 1:n_ps){
     each_tree_pred[,i] <- B_train%*%beta_mean[,i]
     points(x, each_tree_pred[,i]+matrix(model_run$BUGSoutput$mean$beta_intercept, nrow = nrow(B_train)),
            pch=20, col = alpha(i,0.3)) # True line
}
# lines(x, y_train_hat_low, col = "blue", lty = 2) # Predicted low
# lines(x, y_train_hat_up, col = "blue", lty = 2) # Predicted high
points(x, y_train_hat,pch=20) # True line

model_run$BUGSoutput$mean$tau_b

