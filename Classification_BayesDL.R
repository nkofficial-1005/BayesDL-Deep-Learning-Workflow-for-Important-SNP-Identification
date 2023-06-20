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

#Binary phenotype
#Import SNP (or feature) data set for 
#Anthocyanin_22
#####
xTrain <- read.csv("XTrain_Antho.csv")
xTrain <- xTrain[,-1]
xTrain <- as.data.frame(xTrain)
xTest <- read.csv("XTest_Antho.csv")
xTest <- xTest[,-1]
xTest <- as.data.frame(xTest)
#####

#Import response (or phenotype) data set for
#Anthocyanin_22
#####
yTrain <- read.csv("ytrain_Antho.csv")
yTrain <- yTrain[,-1]
yTrain <- as.vector(yTrain)
yTrain <- replace(yTrain, yTrain==1, 2)
yTrain <- replace(yTrain, yTrain==0, 1)
yTest <- read.csv("ytest_Antho.csv")
yTest <- yTest[,-1]
yTest <- as.vector(yTest)
yTest <- replace(yTest, yTest==1, 2)
yTest <- replace(yTest, yTest==0, 1)
######

#Categorical phenotype
#Import feature data set for 
#Germ_22
#####
# xTrain <- read.csv("XTrain_Germ.csv")
# xTrain <- xTrain[,-1]
# xTrain <- as.data.frame(xTrain)
# xTest <- read.csv("XTest_Germ.csv")
# xTest <- xTest[,-1]
# xTest <- as.data.frame(xTest)
#####

#Import response data set for
#Germ_22
#####
# yTrain <- read.csv("ytrain_Germ.csv")
# yTrain <- yTrain[,-1]
# yTrain <- as.vector(yTrain)
# yTest <- read.csv("ytest_Germ.csv")
# yTest <- yTest[,-1]
# yTest <- as.vector(yTest)
######
x <- bind_rows(xTrain, xTest)
y <- as.data.frame(c(yTrain, yTest))
colnames(y) <- "y"
data <- bind_cols(x,y)
######

#Bayesian Neural Network (BNN)

#Classification NN function 
######
sm <- stan_model("nn_class.stan")

fit_nn_cat <- function(x_train, y_train, x_test, y_test, H, n_H, data, method = "optimizing", ...) {
  stan_data <- list(
    N = nrow(x_train),
    P = ncol(x_train),
    x = x_train,
    labels = y_train,
    H = H,
    n_H = n_H,
    N_test = length(y_test)
  )
  if(method == "optimizing") {
    optOut <- optimizing(sm, data = stan_data)
    test_char <- paste0("output_test[",1:length(y_test), ",",rep(1:max(y_train), each = length(y_test)),"]") 
    y_test_pred <- matrix(optOut$par[test_char], stan_data$N_test, max(y_train))
    y_test_cat <- apply(y_test_pred, 1, which.max)
    out <- list(y_test_pred = y_test_pred,
                y_test_cat = y_test_cat,
                conf = table(y_test_cat, y_test),
                fit = optOut)
    return(out)
  } else if(method == "sampling") {
    out <- sampling(sm, data = stan_data, ...)
    y_test_pred <- summary(out, pars = "output_test")$summary
    out <- list(y_test_pred = y_test_pred,
                fit = out)
    return(out)
  } 
}
######

#Fit class NN
#Optimizing the model
fit_opt <- fit_nn_cat(xTest, yTest,xTrain, yTrain, 2, 50, data, method = "optimizing")

#Calculate Performance Metrics
cm <- fit_opt$conf
cm
#accuracy <- sum(cm[1], cm[4]) / sum(cm[1:4])
precision <- cm[4] / sum(cm[4], cm[2])
sensitivity <- cm[4] / sum(cm[4], cm[3])
fscore <- (2 * (sensitivity * precision))/(sensitivity + precision)
auc <- auc(as.matrix(yTest),fit_opt$y_test_cat)

#Sampling from the fitted model
fit_nuts <- fit_nn_cat(xTrain, yTrain, xTest, yTest, 2, 50, method = "sampling", 
                       chains = 4, cores = 4, iter = 2000, warmup=1000)

#Save the fitted model for future use because running takes a while
saveRDS(fit_nuts, 'stan_fit_antho.rds')
fit <-readRDS('stan_fit_antho.rds')

#Find sampler for each chain
sample <- get_sampler_params(fit_nuts$fit, inc_warmup = TRUE)
lapply(sample, summary, digits = 2)
sapply(sample, FUN = colMeans)

#Parameter names of the draws
list_of_draws <- rstan::extract(fit_nuts$fit)
print(names(list_of_draws))

#Check Rhat values
print(fit_nuts$fit)

#Feature selection using BNN
######
#Extract weights associated with each predictor (or SNP)
P = ncol(xTrain)
N = nrow(xTrain)
n_H = 50
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

# Sort the variable importance measures in ascending order
var_imp <- sort(var_imp, decreasing = FALSE)

# Select the top 10 predictors
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
ppc_dens_overlay((x), t(as.matrix(transformed_params$`data_to_hidden_weights[89,5]`)))
ppc_dens_overlay((x), t(as.matrix(transformed_params$`data_to_hidden_weights[15,23]`)))
