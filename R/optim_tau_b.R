library(boot)
library(mvnfast)

data("motor")
x <- as.matrix(motor$times)
y <- as.matrix(motor$accel)
set.seed(42)
n_ <- 100
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_))
colnames(x) <- "x"
colnames(x_new) <- "x"
y <- sin(3*x) + rnorm(n = n_,sd = 0.1)
# y[x<0] <- y[x<0] + 2
# y[x>0] <- y[x>0] - 2


nIknots <- 10
knots <- quantile(as.vector(x), seq(0, 1, length.out = nIknots + 2))[-c(1, nIknots + 2)]
B <- cbind(1, splines::ns(x = x, intercept = FALSE, knots = knots))


nll <- function(dat, x, par,B) {
     tau <- par[1]
     tau_b <- par[2]
     y <- dat
     n <- length(y)
     tryCatch(-mvnfast::dmvn(t(y), rep(0, n), diag(tau, n) + tau_b * tcrossprod(B), log = TRUE),
     error=function(e) -Inf)
     return(-mvnfast::dmvn(t(y), rep(0,n), diag(tau^-1, n) + (tau_b^-1) * tcrossprod(B) , log = TRUE))
}

nll_b0 <- function(dat, x, par) {
     tau <- par[1]
     tau_b0 <- par[2]
     tau_b <- par[3]
     y <- dat
     n <- length(y)
     B_new <- sweep(B, 2, c(tau_b0, rep(tau_b, ncol(B) - 1)), FUN="*", check.margin=FALSE)
     tryCatch(-mvnfast::dmvn(t(y), rep(0, n), diag(tau, n) + tcrossprod(B_new), log = TRUE),
              error=function(e) -Inf)
}

fit <- optim(par = rep(1, 2), fn = nll, dat=y, x = x,
             method = "L-BFGS-B",
             hessian = TRUE,
             lower = rep(0.01, 2))

fit_b0 <- optim(par = rep(1, 3), fn = nll_b0, dat=y, x=x,
                method = "L-BFGS-B",
                hessian = TRUE,
                lower = rep(0.01, 3))
