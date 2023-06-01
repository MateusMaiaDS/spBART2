rm(list=ls())
library(tidyverse)

n_ <- 100
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_*100))
colnames(x) <- "x"
colnames(x_new) <- "x"
# x <- as.data.frame(x)
# x_test <- as.data.frame(x_new)
y <- sin(2*x) + rnorm(n = n_,sd = 0.1)
y[x<0] <- y[x<0] + 2
y[x>0] <- y[x>0] - 2


# plot(x,motrbartmod$y_hat %>% colMeans())
