# BayesDL-Deep-Learning-Workflow-for-Important-SNP-Identification

##  Table of content

- [Abstract](#abstract)
  * [Glossary](#glossary)

- [Prerequisites](#prerequisites)
  * [Data](#data)
  
- [Methodology](#methodology)

- [Instructions to Run Code](#instructions-to-run-code)
  
- [License](#license)
  
- [Contact](#contact)

##  Abstract

Genome-Wide Association Study (GWAS) is the discovery of  an association between certain variations in the genetic code (genome) and a certain physical trait (phenotype). Single Nucleotide Polymorphisms (SNPs) are the most abundant form of common simple DNA variants. In bioinformatics studies, one of the most challenging processes to carry out association tests is finding significant SNPs in high-dimensional data. This problem can be potentially solved by feature selection using statistical and machine learning algorithms. A Bayesian Neural Network (BNN) based workflow, BayesDL, is proposed to identify important SNPs. BayesDL is a cascaded classifier and regressor of an Artificial neural network using Bayesian inference. The model is fitted using R and Stan. Firstly, the whole-genome SNP data undergoes preliminary feature (SNP) selection using hypothesis testing procedures to reduce the size of the data. It is needed to make the data compatible with neural network (NN) architectures (especially BNN which also includes probabilistic distributions) as they require extensive computational resources (more than 32GB RAM when using all SNPs as input layer). Secondly, the BNN model is defined using Stan; the defined model is further used for SNP identification and test performance evaluation. Finally, the existing CNN model is utilized for comparison based on test performance metrics which aids in determining the superiority of BayesDL. The final output (important SNPs) from BayesDL is used for further biological analysis. The following paragraphs provide steps followed to perform preliminary feature selection and define the research work required to train and validate the NN-based models.

BayesDL, used for SNP identification, increases confidence in the selected SNPs by leveraging Bayesian methods and deep learning advantages. Also, by providing evidence against the hypothesis that apriori SNPs are insignificant, the proposed workflow aids in increasing confidence in the selected set of SNPs. Unlike a CNN model, which uses a single value for the weights, BayesDL accounts for the uncertainty of those weights by using probabilistic distribution, thereby improving the confidence in the chosen SNP set. Additionally, using BayesDL is beneficial when working with small sample sizes, as it helps to reduce over-fitting, a common issue evident in CNN. Using the BNN model also results in less number of False negatives as compared to that when using CNN.

Key Words: Genomic Wide Association Study  ·  Single Nucleotide Polymorphism  ·  Feature Selection  ·  Deep Learning  ·  High Dimensional Data.

### Glossary

* SNP - Single Nucleotide Polymorphisms
* GWAS - Genome-Wide Association Studies
* TASSEL - Trait Analysis by aSSociation, Evolution and Linkage
* GAPIT -  Genome Association and Prediction Integrated Tool
* DL - Deep Learning
* NN - Neural Network
* CNN - Convolutional Neural Network
* ANN - Artificial Neural Network
* CoV - Coefficient of Variation
* MCMC - Markov Chain Monte Carlo
* GLM - Generalized Linear Model
  
## Prerequisites

Software: R Version 4.2.2 and R Version 3.6.3

Operating Systems: Linux 5.4.0-135-generic x86_64 and Linux 5.4.0-150-generic x86_64

###  Data

Arabidopsis thaliana  data, AtPolyDB, is used for this study. It is obtained from the easygwas website: https://easygwas.ethz.ch/data/public/dataset/view/1/. The AtPolyDB dataset has 1307 samples with 214051 SNPs (or features) and the F1 data set has 372 samples with 204753 SNPs. The data set contains three files: (a) PED file, (b) PHENO file, and (c) MAP file. The chosen phenotypes had three different data types: (a) Binary (Anthocyanin), (b) Continuous (Width), and (c) Categorical (Germination Days).

##  Methodology

The new workflow includes below steps,

<b>Input:</b>  Genotype .ped file and Phenotype .pheno file  
1. Pre-process the data and perform feature selection using chi-square and ANOVA tests for categorical and continuous phenotypes respectively.
2. Split the data into training (50\%) and testing (50\%) folds.
3. Develop and specify the Stan model as follows: 
a. The RStan model is developed to form data, transformed data, parameters, function, transformed parameters, model, and generated quantities block. Save the files as nn_reg.stan and nn_class.stan for further analysis.  
b. Specify the stan models in R using nn_reg.stan and nn_class.stan files. Define a function to compile the stan model and get the desired output as follows:
	1. Define stan data in the function.  
	2. Call Stan's optimizing methods to obtain point estimates by optimizing the posterior for the model. 
	3. Call a Stan's NUTS sampler to draw posterior samples from the model.
4. Fit train and test sets (and vice-versa) in the function for training the model and making predictions using $2$ hidden layers with $50$ neurons in each layer respectively. 
5. Record the two test performance metrics for each phenotype from the optimization method.
6. Rank the input SNPs by generating the samples of weights (weights ranked in increasing value of CoV) corresponding to predictors using the sampling function.

<b>Output:</b>  The top $10$ significant SNPs for each phenotype.

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Contact
You can ask questions to [Nikita Kohli](mailto:nikita.datascience@gmail.com).
