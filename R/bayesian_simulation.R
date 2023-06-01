# Getting a simulation over an x
bayes_sim <- function(n_, nIknots_,seed_,
                      tau_b,tau,n_post_,n_tree_ = 1)
{
set.seed(seed_)
x <- runif(n = n_)
knots <- quantile(x,seq(0,1,length.out = nIknots_+2))[-c(1,nIknots_+2)]
B <- splines::ns(x = x,knots = knots,intercept = FALSE)
B <- cbind(1,B) # Adding the intercept column
kappa <- 2
tau_b0 <- 4*(kappa^2)*n_tree_
# Generating samples from the model
n_post_ <- 2000
beta_post_ <- matrix(NA,nrow = n_post_,ncol = ncol(B))
y <- matrix(NA,nrow = n_post_,ncol = n_)
beta_cov_ <- diag(tau_b^(-1),nrow = ncol(B))
beta_cov_[1,1] <- tau_b0^(-1)

for(i in 1:n_post_){
     beta_post_[i,] <- mvnfast::rmvn(n = 1,mu = rep(0,ncol(B)),sigma = beta_cov_)
     y[i,] <- tcrossprod(beta_post_[i,],B) + rnorm(n = n_,mean = 0,sd = sqrt(tau^-1))
}

return(list(x = x,
            y = y))

}

