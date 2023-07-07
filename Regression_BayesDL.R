# Clear environment
##############
rm(list=ls())
gc()
options(warn=-1)

#####
library(rstan)
library(magrittr)
library(caret)
library(tidyverse)
library(dplyr)
library(pROC)
library(ROCR)
library(bayesplot)
library(shinystan)
######

#Continuous phenotype
#Import SNP (or feature) data set for 
#Width_22
######
xTrain <- read.csv("XTrain_Width.csv")
xTrain <- xTrain[,-1]
xTrain <- as.data.frame(xTrain)
xTest <- read.csv("XTest_Width.csv")
xTest <- xTest[,-1]
xTest <- as.data.frame(xTest)
######

#Import response (or phenotype) data set for
#Width_22
######
yTrain <- read.csv("ytrain_Width.csv")
yTrain <- yTrain[,-1]
yTrain <- as.vector(yTrain)
yTest <- read.csv("ytest_Width.csv")
yTest <- yTest[,-1]
yTest <- as.vector(yTest)
######
x <- bind_rows(xTrain, xTest)
y <- as.data.frame(c(yTrain, yTest))
colnames(y) <- "y"
data <- bind_cols(x,y)
######

#Bayesian Neural Network (BNN)

#Regression NN function
###### 
sm_reg <- stan_model("nn_reg.stan")

fit_nn_reg <- function(x_train, y_train, x_test, y_test, H, n_H, data, method = "optimize", ...) {
  stan_data <- list(
    N = nrow(x_train),
    P = ncol(x_train),
    x = x_train,
    y = y_train,
    H = H,
    n_H = n_H,
    N_test = length(y_test),
    y_test = y_test
  )
  if(method == "optimize") {
    optOut <- optimizing(sm_reg, data = stan_data)
    test_char <- paste0("output_test[", 1:stan_data$N_test, "]")
    y_test_pred <- optOut$par[test_char]
    mse <- mean((y_test_pred - y_test)^2)
    correlation <- cor(y_test_pred,y_test)
    rsq <- (correlation^2)
    out <- list(y_test_pred = y_test_pred,
                sigma = optOut$par["sigma"],
                mse  = mse,
                rsq = rsq,
                fit = optOut)
    return(out)
  } else {
    if(method == "sampling") {
      out <- sampling(sm_reg, data = stan_data, ...)
    } else if (method == "vb") {
      out <- vb(sm_reg, data = stan_data, pars = c("output_test", "sigma", "output_test_rng"), ...)
    }
    y_test_pred <- summary(out, pars = "output_test")$summary
    sigma <- summary(out, pars = "sigma")$summary
    out <- list(y_test_pred = y_test_pred,
                sigma = sigma,
                fit = out)
    return(out)
  }
}
######

#Fit regression NN
#Optimizing the model
fit_opt <- fit_nn_reg(xTest, yTest, xTrain, yTrain, 2, 50, data, method = "optimize")

#Calculate Performance Metrics
#RMSE <- sqrt(fit_opt$mse)
#R2 <- fit_opt$rsq

#Sampling from the fitted model
fit_nuts <- fit_nn_reg(xTrain, yTrain, xTest, yTest, 2, 50, method = "sampling",
                       chains = 4, cores = 4, iter = 2000, warmup=1000)

#Save the fitted model for future use because running takes a while
saveRDS(fit_nuts, 'stan_fit_width.rds')
fit <-readRDS('stan_fit_width.rds')

#Find sampler for each chain
sample <- get_sampler_params(fit_nuts$fit, inc_warmup = TRUE)
lapply(sample, summary, digits = 2)
sapply(sample, FUN = colMeans)

#Parameter names of the draws
list_of_draws <- rstan::extract(fit_nuts$fit)
print(names(list_of_draws))

#Check the Rhat values
print(fit_nuts$fit)

#Feature selection using BNN
######
#Extract weights associated with each predictor
P = ncol(xTrain)
N = nrow(xTrain)
n_H = 50 #Number of nodes in each hidden layer
wt_samples <- matrix(NA, nrow = P, ncol = n_H)

# Define the names of the parameters which are to be extracted
param_names <- c()
for (i in 1:P){
  for(j in 1:n_H){
    wt_name <- paste0("data_to_hidden_weights[", i, ",", j, "]")
    param_names <- c(param_names, wt_name)
  }
}

# Extract the parameter samples and store them in a matrix
wt_samples <- matrix(NA, nrow = 4000, ncol = length(param_names))
col_name <- list()
for (i in 1:length(param_names)) {
  wt_name <- param_names[i]
  wt_samples[, i] <- rstan::extract(fit_nuts$fit, pars = wt_name)[[1]]
  wt_samples <- as.data.frame(wt_samples)
  col_name = append(col_name,param_names[i])
}
colnames(wt_samples) <- col_name

# Compute the posterior means and standard deviations for each predictor (or SNP)
wt_means <- colMeans(wt_samples)
wt_sds <- apply(wt_samples, 2, sd)

#Remove weight means that are equal to zero
lapply(wt_means, function(x) {x[x!=0]})

# Compute the variable importance measures
var_imp <- wt_sds/ abs(wt_means)

# Sort the variable importance measures in descending order
var_imp <- sort(var_imp, decreasing = FALSE)

# Select the top 10 SNPs
top_vars <- names(var_imp)[1:10]
df <- as.data.frame(top_vars)

# Extract the first number before the comma from each row
df$first_num <- str_extract(as.character(df$top_vars), "\\d+")
df$top_vars
num <- as.integer(df$first_num)
selectedSNPs <- list()
for (i in num) {
  selectedSNPs <- append(selectedSNPs, colnames(xTrain)[i])
}
selectedSNPs <- as.data.frame(selectedSNPs)
######

#MCMC Diagnosis for top 2 SNPs
#######
#Nuts and posterior parameters
np_draws <- nuts_params(fit_nuts$fit)
posterior_draws <- as.array(fit_nuts$fit)

# Pairs Plot
color_scheme_set("red")
mcmc_pairs(posterior_draws, np = np_draws, pars = top_vars[1:2],
           off_diag_args = list(size = 0.75),
           main = top_vars[1:2])

#Trace Plot
color_scheme_set("green")
mcmc_trace(posterior_draws, pars = top_vars[1:2], np = np_draws) +
  xlab("Post-warmup iteration")

#Acf plot
color_scheme_set("mix-brightblue-gray")
mcmc_acf(posterior_draws, pars = top_vars[1:2], lags = 35)

#Posterior Predictive checks
params <- rstan::extract(fit_nuts$fit, pars = top_vars[1:2])

#Prior Distribution
x <- rnorm(4000,mean = 0, sd = 1)
x <- as.data.frame(x)
x <- as.numeric(unlist(x))

#PPC Plot
color_scheme_set("viridis")
par(mfrow=c(1,2))
ppc_dens_overlay((x), t(as.matrix(transformed_params$`data_to_hidden_weights[13,25]`)))
ppc_dens_overlay((x), t(as.matrix(transformed_params$`data_to_hidden_weights[19,4]`)))
