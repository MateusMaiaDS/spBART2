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
n_ <- 500 # Number of observations

# Simulation 2

no_interaction_friedman <- function (n, sd = 1)
{
     x <- matrix(runif(n), ncol = 1)
     y <- 10 * sin(pi * x[, 1] )
     # y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
     if (sd > 0) {
          y <- y + rnorm(n, sd = sd)
     }
     list(x = x, y = y)
}

# Simulation 1
sd_ <- 5
fried_sim <- no_interaction_friedman(n = n_,sd = sd_)
x <- fried_sim$x
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
                             nrow(knots)+3, # +3 here because is a natural spline
                             ncol(x_train[,continuous_vars, drop = FALSE])))

B_test_arr <- array(data = NA,
                    dim = c(nrow(x_test),
                            nrow(knots)+3,  # +3 here because is a natural spline
                            ncol(x_test[,continuous_vars, drop = FALSE])))


# Getting min and max
min_x <- x_min <- apply(as.matrix(x_train),2,min)
max_x <- x_max <- apply(as.matrix(x_train),2,max)


# Getting the internal knots
nIknots <- 50
dif_order <- 2
knots <- apply(x_train,
               2,
               function(x){quantile(x,seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]})



B_train_arr <- array(data = NA,
                     dim = c(nrow(x_train),
                             nrow(knots)+3, # +3 here because is a natural spline
                             ncol(x_train[,continuous_vars, drop = FALSE])))

B_test_arr <- array(data = NA,
                    dim = c(nrow(x_test),
                            nrow(knots)+3,  # +3 here because is a natural spline
                            ncol(x_test[,continuous_vars, drop = FALSE])))

Z_train_arr <- array(data = NA,
                     dim = c(nrow(x_train),
                             nrow(knots)+3,
                             ncol(x_train[,continuous_vars, drop = FALSE]))) # correcting the new dimension by P

Z_test_arr <- array(data = NA,
                    dim = c(nrow(x_test),
                            nrow(knots)+3,
                            ncol(x_test[,continuous_vars, drop = FALSE])))# correcting the new dimension by P


# Creating the natural B-spline for each predictor
for(i in 1:length(continuous_vars)){
     B_train_obj <- splines::bs(x = x_train[,continuous_vars[i]],knots = knots[,continuous_vars[i]],
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
                                  nrow(knots)+3-dif_order,
                                  ncol(x_train[,continuous_vars, drop = FALSE]))) # correcting the new dimension by P

     Z_test_arr <- array(data = NA,
                         dim = c(nrow(x_test),
                                 nrow(knots)+3-dif_order,
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


# Setting other parameters
nsigma <- naive_sigma(x = x_train,y = y_train)
df <- 3
# Calculating tau hyperparam
a_tau <- df/2
sigquant <- 0.9

tau_b_0 <- 4*(2^2)*1/(diff(range(y_train))^2)

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

     mean_y[t] = inprod(B_one[t,1:N_knots],beta_one[1:N_knots]) #+ inprod(B_two[t,1:N_knots],beta_two[1:N_knots]) + inprod(B_three[t,1:N_knots],beta_three[1:N_knots]) + inprod(B_four[t,1:N_knots],beta_four[1:N_knots]) #+ inprod(B_five[t,1:N_knots],beta_five[1:N_knots])

     y[t] ~ dnorm(sum(mean_y[t]) + gamma, tau)
   }



    # RW prior on beta
    gamma ~ dnorm(0,tau_b_0)

    for(i in 1:N_knots){
               beta_one[i] ~ dnorm(0,tau_b[1])
               # beta_two[i] ~ dnorm(0,tau_b[2])
               # beta_three[i] ~ dnorm(0,tau_b[3])
               # beta_four[i] ~ dnorm(0,tau_b[4])
               # beta_five[i] ~ dnorm(0,tau_b[5])
    }

    # Priors on beta values
    tau ~ dgamma(a_tau, d_tau)
    for(k in 1:1){
         tau_b[k] ~ dgamma(0.5 * nu, 0.5 * delta[k] * nu)
         delta[k] ~ dgamma(a_delta, d_delta)
    }
  }
  "

# for(i in 1:length(nu_range)){

# Set up the data
model_data <- list(N = nrow(Z_train_arr[,,1]),
                   y = as.vector(y_train),
                   B_one = Z_train_arr[,,1],
                   # B_two = Z_train_arr[,,2],
                   # B_three = Z_train_arr[,,3],
                   # B_four = Z_train_arr[,,4],
                   # B_five = Z_train_arr[,,5],
                   N_knots = ncol(Z_train_arr[,,1]),
                   # diag_ = diag(nrow = ncol(B_train_arr[,,1])),
                   a_tau = a_tau,
                   d_tau = d_tau,
                   a_delta = 0.0001,
                   d_delta = 0.0001,
                   tau_b_0 = tau_b_0,
                   nu  = 2)

# Choose the parameters to watch
model_parameters <- c("beta_one", "tau", "tau_b","delta","gamma")

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
beta_post_one <- model_run$BUGSoutput$sims.list$beta_one
beta_quantile_one <- apply(beta_post_one, 2, quantile, prob = c(0.25, 0.5, 0.75))
# beta_post_two <- model_run$BUGSoutput$sims.list$beta_two
# beta_quantile_two <- apply(beta_post_two, 2, quantile, prob = c(0.25, 0.5, 0.75))
# beta_post_three <- model_run$BUGSoutput$sims.list$beta_three
# beta_quantile_three <- apply(beta_post_three, 2, quantile, prob = c(0.25, 0.5, 0.75))
# beta_post_four <- model_run$BUGSoutput$sims.list$beta_four
# beta_quantile_four <- apply(beta_post_four, 2, quantile, prob = c(0.25, 0.5, 0.75))
# beta_post_five <- model_run$BUGSoutput$sims.list$beta_five
# beta_quantile_five <- apply(beta_post_five, 2, quantile, prob = c(0.25, 0.5, 0.75))
gamma_post <- model_run$BUGSoutput$sims.list$gamma
gamma_quantile <- apply(gamma_post, 2, quantile, prob = c(0.25, 0.5, 0.75))


# New prediction
# Z_test <- predict(Z_train,x_test)
y_train_hat <- Z_train_arr[,,1]%*%beta_quantile_one[2,] + gamma_quantile[2,]#+ Z_train_arr[,,2]%*%beta_quantile_two[2,]+
     # Z_train_arr[,,3]%*%beta_quantile_three[2,] + Z_train_arr[,,4]%*%beta_quantile_four[2,]+
     # Z_train_arr[,,5]%*%beta_quantile_five[2,] +
     # gamma_quantile[2,]


# Running BART
bartmod <- dbarts::bart(x.train = x_train,y.train = unlist(c(y_train)),ntree = 200,x.test = x_test,keeptrees = TRUE)

# Convergence plots
par(mfrow = c(1,2))
plot(model_run$BUGSoutput$sims.list$tau,type = "l", main = expression(tau),ylab=  "")
plot(bartmod$sigma^-2, type = "l", main = paste0("BART: ",expression(tau)),ylab=  "")

par(mfrow = c(1,2))
plot(y_train_hat,y, main = 'mpsBART', xlab = "mpsBART pred", ylab = "y")
plot(bartmod$yhat.train.mean,y, main = "BART", xlab = "BART pred", ylab = "y")


par(mfrow = c(1,1))
plot(c(unlist(x_train)),y)
points(c(unlist(x_train)),y_train_hat,pch=20)
points(c(unlist(x_train)),10*sin(pi*c(unlist(x_train))), col = alpha("orange",0.5),pch = 20)

# Generating new samples
new_sample <- no_interaction_friedman(n = 100,sd = sd_)

# Calculations in the splines model
B_test_ <- predict(B_train_obj,new_sample$x)
if(dif_order!=0){
     Z_test_ <- B_test_%*%crossprod(D,solve(tcrossprod(D)))
     y_test_hat <- Z_test_%*%beta_quantile_one[2,] + gamma_quantile[2,]#+ Z_train_arr[,,2]%*%beta_quantile_two[2,]+

} else{
     y_test_hat <- B_test_%*%beta_quantile_one[2,] + gamma_quantile[2,]#+ Z_train_arr[,,2]%*%beta_quantile_two[2,]+

}

# Calculations in the bart mod
y_test_bart <- colMeans(predict(bartmod,new_sample$x))

plot(new_sample$x,y_test_hat)
rmse(x = y_test_hat,y = new_sample$y)


plot(new_sample$x,y_test_bart)
rmse(x = y_test_bart,y = new_sample$y)
