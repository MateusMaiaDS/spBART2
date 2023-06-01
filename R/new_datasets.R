# nIknots <- 10
# knots <- quantile(unlist(c(x)),seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]
# B <- splines::ns(x = as.matrix(x),intercept = FALSE,knots = knots)
#
# # Faithfull dataset
# data("faithful")
# x <- faithful$waiting %>% as.matrix()
# y <- faithful$eruptions %>% as.matrix()
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
# x_new <- x
#
#
# # cars
# x <- cars$speed %>% as.matrix()
# y <- cars$dist %>% as.matrix()
# colnames(x) <- "x"
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
# x_new <- x
#
#
# # mtcars
# x <- mtcars$mpg %>% as.matrix()
# y <- mtcars$hp %>% as.matrix()
# colnames(x) <- "x"
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
# x_new <- x
#
#
# # pressure
# x <- pressure$temperature %>% as.matrix()
# y <- pressure$pressure %>% as.matrix()
# colnames(x) <- "x"
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
# x_new <- x


# Diabetes
# diabetes <- read_csv("data/diabetes/diabetes.csv")
# x <- diabetes$age %>% as.matrix()
# y <- diabetes$cpep
#
# x <- diabetes$base %>% as.matrix()
# y <- diabetes$cpep
#
#
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
#
#
# # Fetal growth
# fetal_growth <- read_csv("data/fetal_growth/fetal_growth.csv")
# x <- fetal_growth$gawks %>% as.matrix()
# y <- fetal_growth$hemi
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
# # Nerve
# nerve <- read_csv("data/nerve/nerve.csv")
# x <- nerve$vel %>% as.matrix()
# y  <- nerve$age
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
# # body fat
# res_bodyfat <- read_csv("data/res_bodyfat/res_bodyfat.csv")
# x <- res_bodyfat$bmi %>% as.matrix()
# y <- res_bodyfat$pbfm
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
# # tricpes
# triceps <- read_csv("data/triceps/triceps.csv")
# x <- triceps$age %>% as.matrix()
# y <- triceps$lntriceps
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
# # cps71
# library(crs)
# data("cps71")
# x <- cps71$logwage %>% as.matrix()
# y <- cps71$age
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
# # Engel
# library(quantreg)
# data("engel")
# x <- engel$income %>% as.matrix()
# y <- engel$foodexp
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
# # Corbet
# library(VGAM)
# data("corbet")
# x <- corbet$ofreq %>% as.matrix()
# y <- corbet$species
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
#
# # Hormone
# library(VGAM)
# data("hormone")
# x <- hormone$X %>% as.matrix()
# y <- hormone$Y
# x_new <- seq(min(x),max(x),length.out = 100) %>% as.matrix()
# colnames(x) <- colnames(x_new) <- "x"
#
# # Motor
# library(boot)
# data("motor")
# x <- motor$times %>% as.matrix
# y <- motor$accel %>% as.matrix()
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
#
# # # Faithfull dataset
# data("faithful")
# x <- faithful$waiting %>% as.matrix()
# y <- faithful$eruptions %>% as.matrix()
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
#
# # cars
# x <- cars$speed %>% as.matrix()
# y <- cars$dist %>% as.matrix()
# colnames(x) <- "x"
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
#
# # # mtcars
# x <- mtcars$mpg %>% as.matrix()
# y <- mtcars$hp %>% as.matrix()
# colnames(x) <- "x"
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"
#
# #Abalone
# abalone <- read.table(file = "data/abalone.data",sep = ",")
# # abalone_names <- read.table(file = "data/abalone/abalone.names")
# x <- abalone$V2 %>% as.matrix()
# y <- abalone$V8
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"

# Simulation review
# set.seed(42)
# n_ <- 400
# x <- (0:(n_-1) / (n_-1)) %>% as.matrix()
# f <- -3.5+0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
# y <- f + rnorm(n_, 0, sd = 2)
# x_new <- seq(min(x),max(x),length.out = 1000) %>% as.matrix()
# colnames(x) <- "x"
# colnames(x_new) <- "x"

# Getting the simulation
# data_ <- bayes_sim(n_ = 2000,nIknots_ = 10,
#           seed_ = 42,tau_b = 0.01,
#           tau = 1,n_post_ = 2000,n_tree_ = 1)
# x <- data_$x
# y <- colMeans(data_$y)
