################################################################################
######                                                                   #######
######                    Code for generating example data               #######
######                                                                   #######
################################################################################

source("/Users/mike/Documents/Mike/Projects/Reduced Rank/Simulation/new_code/Gen_data.R")
library(MASS)


rho <- 0.1; b0 <- 0.2;
n <- 50
p <- q <- m <- 20
r <- 3

X = getX(rho, n, p, q, high = F)
A = getA(b0, p, r, m)
Y = X %*% A + getE(n, m)

usethis::use_data(X)
usethis::use_data(Y)
