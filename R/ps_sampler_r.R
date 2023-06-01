rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/spbart.cpp")
source("R/other_functions.R")
source("R/wrap_bart.R")
source("R/bayesian_simulation.R")
n_ <- 500
set.seed(42)

# Simulation 1
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_))
colnames(x) <- "x"
colnames(x_new) <- "x"
y <- sin(3*x) + rnorm(n = n_,sd = 0.1)

# Creating the basis matrix
nIknots <- 50
knots <- quantile(x[,1],seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]
B_train <- splines::bs(x = x,knots = knots)
B_test <- predict(B_train,x_new)
p <- ncol(B_train)
D <- matrix(0, p-2, p)
# Populating D
for(j in 1:(p-2)) {
     D[j,j:(j+2)] <- c(1,-2,1)
}
P <- crossprod(D)
# P <- diag(nrow = ncol(B_train))

# Creating each parameter sampler
betas_sampler <- function(B_train,P,y,tau_b,tau_b_zero,tau){

     n <- (length(y))
     ones <- matrix(1,nrow = n)

     # Creating the elements
     Gamma_inv <- solve(crossprod(B_train) + (tau_b/tau)*P - (1/(n+tau_b_zero/tau))*crossprod(B_train,ones)%*%crossprod(ones,B_train))

     beta_mean <- Gamma_inv%*%(crossprod(B_train,y))
     beta_var <- (1/tau)*Gamma_inv

     beta_sample <- mvnfast::rmvn(n = 1,mu = beta_mean,sigma = beta_var)

     return(beta_sample)
}

beta_zero_sampler <- function(B_train,P,betas,tau_b_zero,tau){

     n <- (length(y))
     ones <- matrix(1,nrow = n)
     betas <- matrix(betas,nrow = ncol(B_train))
     beta_zero_mean <- (1/(n+tau_b_zero/tau))*(sum(y)-crossprod(betas,crossprod(B_train,ones)))
     beta_zero_mean_var <- (1/(n+tau_b_zero/tau))*(1/tau)

     beta_zero_sample <- rnorm(n = 1,mean = beta_zero_mean,sd = sqrt(beta_zero_mean_var))
}

tau_b_sampler <- function(betas,nu,delta){

     tau_b_shape <- 0.5*length(betas)+0.5*nu
     tau_b_rate <- 0.5*(betas%*%tcrossprod(P,betas))+0.5*delta*nu

     tau_b_sample <- rgamma(n = 1,shape = tau_b_shape,rate = tau_b_rate)
     return(tau_b_sample)
}

delta_sampler <- function(nu,tau_b,a_delta,d_delta){

     delta_shape <- 0.5*nu + a_delta
     delta_rate <- 0.5*nu*tau_b + d_delta
     delta_sample <- rgamma(n = 1,shape = delta_shape,rate = delta_rate)
     return(delta_sample)
}

tau_sampler <- function(B_train,betas,y,a_tau,d_tau){
     y_hat <- tcrossprod(B_train,betas)
     return(rgamma(n = 1,shape = 0.5*length(y)+a_tau,rate = 0.5*crossprod(y)+d_tau))
}

# Setting all hyperparameters and initial values
x_min <- min(x); y_min <- min(y)
x_max <- max(x); y_max <- max(y)
x <- normalize_covariates_bart(y = x,a = x_min,b = x_max)
y <- normalize_bart(y = y,a = y_min,b = y_max)
tau_b <- 16
tau <- 1
tau_b_zero <- 4*(2^2)*1
a_delta <- d_delta <- 0.0001
delta <- 0.0001
nu <- 2
df <- 3
# Calculating tau hyperparam
a_tau <- df/2
sigquant <- 0.9
nsigma <- naive_sigma(x = x,y = y)
# Calculating lambda
qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
lambda <- (nsigma*nsigma*qchi)/df
d_tau <- (lambda*df)/2

# Setting MCMC parameters
n_mcmc <- 2500
n_burn <- 500
n_post <- n_mcmc-n_burn

# ================
# INITIALISING THE MODEL
# ================
beta_post <- matrix(data = NA,nrow = n_post,ncol = ncol(B_train))
y_hat_post <- matrix(data = NA,nrow = n_post,ncol = length(y))

beta_zero_post <- numeric()
tau_post <- numeric()
tau_b_post <- numeric()
delta_post <- numeric()
curr <- 1

for(i in 1:n_mcmc){

     betas <- betas_sampler(B_train = B_train,P = P,y = y,tau_b = tau_b,tau_b_zero = tau_b_zero,tau = tau)
     beta_zero <- beta_zero_sampler(B_train = B_train,P = P,betas = betas,tau_b_zero = tau_b_zero,tau = tau)
     tau_b <- tau_b_sampler(betas = betas,nu = nu,delta = delta)
     delta <- delta_sampler(nu = nu,tau_b = tau_b,a_delta = a_delta,d_delta = d_delta)
     tau <- tau_sampler(y = y,a_tau = a_tau,d_tau = d_tau,B_train = B_train,betas = betas)

     if(i > n_burn){

          beta_post[curr,] <- betas
          beta_zero_post[curr] <- beta_zero
          tau_b_post[curr] <- tau_b
          delta_post[curr] <- delta
          tau_post[curr] <- tau
          y_hat_post[curr,] <- tcrossprod(B_train,betas) + beta_zero
          curr = curr+1



     }

}

y_hat_post <- unnormalize_bart(z = y_hat_post,a = y_min,b = y_max)


# Plotting the predictions
par(mfrow=c(1,1))
plot(x,unnormalize_bart(y,a = y_min,b = y_max),main = "P-Splines robust priors")
quantiles_y_hat <- apply(y_hat_post,2,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
lines(x,colMeans(y_hat_post),col = "blue")
lines(x,quantiles_y_hat[1,],lty = "dashed", col = "blue")
lines(x,quantiles_y_hat[3,],lty = "dashed", col = "blue")

# Traceplots
par(mfrow=c(2,2))
plot(beta_zero_post,type = "l", main = expression(beta[0]))
plot(tau_b_post,type = "l", main = expression(tau[b]))
plot(delta_post,type = "l", main = expression(delta))
plot(tau_post,type = "l", main = expression(tau))
