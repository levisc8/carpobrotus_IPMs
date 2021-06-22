# Dependencies

library(fs) # file manipulation
library(ggplot2) # plotting 
library(dplyr) # data manipulation
library(tidyr)
library(gridExtra) # plotting
library(rlang)
library(mgcv)
library(Rcpp)
library(ipmr)
library(nlme)


get_gg_legend<-function(plot){
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


sourceCpp(file = 'Ana_Israel_IPM/Cpp/cpp_utils.cpp')

# Top level wrappers for sensitivity and elasticity computations.

sensitivity <- function(K, h, level = c('kernel', 'vr', 'param'),
                        ...) {
  
  switch(level,
         "kernel" = .kernel_sens(K, h),
         "vr"     = .vr_sens(K, h, ...))
}

elasticity <- function(K, h, level = c('kernel', 'vr', 'param'),
                       ...) {
  
  switch(level,
         "kernel" = .kernel_elas(K, h),
         "vr"     = .vr_elas(K, h, ...))
}


# Sensitivity for kernel, vital rates, and parameters
.kernel_sens <- function(K, h) {
  w <- Re(eigen(K)$vectors[ , 1])
  v <- Re(eigen(t(K))$vectors[ , 1])
  
  out <- outer(v, w) / sum(v * w * h)
  
  return(out)
}


# Elasticity for kernel, vital rates, and parameters

.kernel_elas <- function(K, h) {
  
  sens <- sensitivity(K, h, 'kernel')
  
  lambda <- Re(eigen(K)$values[1])
  
  out <- sens * (K / h) / lambda
  
  return(out)
}

