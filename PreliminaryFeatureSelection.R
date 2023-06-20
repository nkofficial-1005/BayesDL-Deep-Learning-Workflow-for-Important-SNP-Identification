##############
# Clear environment
rm(list=ls())
gc()
options(warn=-1)

##############
# Import libraries
library(dplyr)
library(ROCR)
library(logr)
library(matrixStats)
library(Matrix)
library(glmnet)
library(SGL)
library(pROC)
library(MLmetrics)
library(gglasso)
library(caTools)
library(caret)
library(tidyverse)
library(foreach)
library(doParallel)
library(grpreg)
##############
# Import data
Geno <- read.pedfile("genotype.ped") 
char.pheno <- read.table("phenotypes.pheno", header = TRUE, stringsAsFactors = FALSE, sep = " ")

##############
# Data Pre-processing
set.seed(100)
# Re-code the data in ped file
Geno[Geno == 'A'] <- 0  # Converting A to 0
Geno[Geno == 'T'] <- 1  # Converting T to 1
Geno[Geno == 'G'] <- 2  # Converting G to 2
Geno[Geno == 'C'] <- 3  # Converting C to 3

#Convert the phenotype to matrix
y <- matrix(char.pheno$Anthocyanin_22) #Change the phenotype accordingly
rownames(y) <- char.pheno$IID
index <- !is.na(y)
y <- y[index, 1, drop = FALSE]

#Imputing SNP null values
for (j in 1:ncol(Geno)) {
  Geno[, j] <- ifelse(is.na(Geno[, j]), mean(Geno[, j], na.rm = TRUE), Geno[, 
                                                                            j])
}
Geno_y <- Geno[index, ] 
pheno_final <- data.frame(famid = rownames(y), y = y)

df <- merge(Geno_y, pheno_final, by = 'famid')
df_final <- df[, 7:204760] #comment here what is done
df_final <- sapply(df_final, as.numeric)
df_final <- df_final[sample(nrow(df_final)),]
n <- dim(df_final)[[1]] 
d <- dim(df_final)[[2]] 
##############

#Preliminary Feature Screening
###############
tab = table(df_final[,1], df_final[,d]) 
pvals = rep(NA, d-1)

# Get the chi-squared test for categorical/binary
for(k in 1:d){
  tab = table(df_final[,k], df_final[,d])
  pvals[k] = chisq.test(tab)$p.value
} 

#Anova test for continuous
for(k in 1:d){
  pvals[k] = summary(aov(df_final[,d]~df_final[,k]))[[1]][["Pr(>F)"]][1]
} 

#Finding quantiles
quantile(pvals)

# Identify ones that passed the filter.
a1 <- which(pvals < 0.10)
a2 <- which(pvals < 0.05)
a3 <- which(pvals < 0.01)
a4 <- which(pvals<0.001) #1 in thousand
a5 <- which(pvals<0.0001) #1 in 10 thousand
a6 <- which(pvals<1e-17)

#Create a new data set for filtered SNPs
col_name <- list()
df_new <- matrix(nrow=n,ncol = 0)
df_new <- data.frame(df_new)
for (i in a6) {
  df_new = df_new %>% add_column(df_final[,i])
  col_name = append(col_name,colnames(df_final)[i])
}
colnames(df_new) <- col_name
d1 <- dim(df_new)[[2]]

#Save the filtered data set
write.csv(df_new, "Antho.csv")

#Splitting data in 50% train and 50% test sets
nobs <- nrow(df_new)
id <- createFolds(rowMeans(df_new), k=5, list=F)
training.id <- sample(seq_len(nobs), size = 0.5 * nobs)
testData <- df_new[-training.id, ]
trainData <- df_new[training.id, ]
data.train <- as.data.frame(trainData)
data.test <- as.data.frame(testData)

#Save the train and test sets for further analysis of Neural Networks
x_test <- as.matrix(data.test[,1:(d1-1)])
write.csv(x_test, "XTest_Antho.csv")
y_test <- as.matrix(data.test[, d1])
write.csv(y_test, "ytest_Antho.csv")
x_train <- as.matrix(data.train[,1:(d1-1)])
write.csv(x_train, "XTrain_Antho.csv")
y_train <- as.matrix(data.train[, d1])
write.csv(y_train, "ytrain_Antho.csv")
