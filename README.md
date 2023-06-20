# BayesDL-Deep-Learning-Workflow-for-Important-SNP-Identification

#  Table of content

* [Abstract](#abstract)

* [Technologies](#technologies)

* [Data](#data)

* [Step by step implementation](#step-by-step-implementation)

#  Abstract

Genome-Wide Association Study (GWAS) is the discovery of  an association between certain variations in the genetic code (genome) and a certain physical trait (phenotype). Single Nucleotide Polymorphisms (SNPs) are the most abundant form of common simple DNA variants. In bioinformatics studies, one of the most challenging processes to carry out association tests is finding significant SNPs in high-dimensional data. This problem can be potentially solved by feature selection using statistical and machine learning algorithms. A Bayesian Neural Network (BNN) based workflow, BayesDL, is proposed to identify important SNPs. BayesDL is a cascaded classifier and regressor of an Artificial neural network using Bayesian inference. The model is fitted using R and Stan. Firstly, the whole-genome SNP data undergoes preliminary feature (SNP) selection using hypothesis testing procedures to reduce the size of the data. It is needed to make the data compatible with neural network (NN) architectures (especially BNN which also includes probabilistic distributions) as they require extensive computational resources (more than 32GB RAM when using all SNPs as input layer). Secondly, the BNN model is defined using Stan; the defined model is further used for SNP identification and test performance evaluation. Finally, the existing CNN model is utilized for comparison based on test performance metrics which aids in determining the superiority of BayesDL. The final output (important SNPs) from BayesDL is used for further biological analysis. The following paragraphs provide steps followed to perform preliminary feature selection and define the research work required to train and validate the NN-based models.

BayesDL, used for SNP identification, increases confidence in the selected SNPs by leveraging Bayesian methods and deep learning advantages. Also, by providing evidence against the hypothesis that apriori SNPs are insignificant, the proposed workflow aids in increasing confidence in the selected set of SNPs. Unlike a CNN model, which uses a single value for the weights, BayesDL accounts for the uncertainty of those weights by using probabilistic distribution, thereby improving the confidence in the chosen SNP set. Additionally, using BayesDL is beneficial when working with small sample sizes, as it helps to reduce over-fitting, a common issue evident in CNN. Using the BNN model also results in less number of False negatives as compared to that when using CNN.

Key Words: Genomic Wide Association Study  ·  Single Nucleotide Polymorphism  ·  Feature Selection  ·  Deep Learning  ·  High Dimensional Data.

#  Technologies

Software: R Version 4.2.2 and R Version 3.6.3

Operating Systems: Linux 5.4.0-135-generic x86_64 and Linux 5.4.0-150-generic x86_64

Cloud Servers: TRU Data Science and Compute Canada

#  Data

Two  Arabidopsis thaliana  data, AtPolyDB and F1, are used for this study. They are obtained from easygwas websites: https://easygwas.ethz.ch/data/public/dataset/view/1/ and https://easygwas.ethz.ch/data/public/dataset/view/42/. The AtPolyDB dataset has 1307 samples with 214051 SNPs (or features) and the F1 data set has 372 samples with 204753 SNPs. Both data sets contain three files: (a) PED file, (b) PHENO file, and (c) MAP file. The chosen phenotypes had three different data types: (a) Binary (Anthocyanin), (b) Continuous (Width and DTF), and (c) Categorical (Germination Days).

#  Step-by-step implementation

The new pipeline includes below steps,

<b>Input:</b>  Genotype .ped file and Phenotype .pheno file  
1. Re-code the chromosomal nucleotide to numeric values to form a binary marker data followed by creating a design matrix of dimensions n  ×  p.
2. Remove null values from the phenotype data and match them with marker data.  
3. Impute SNPs with null values with the mean across all samples in the marker data set.  
4. For a 5-fold CV, repeat the steps  
a. Split the data into training (80%) and testing (20%) folds.  
b. Use glmnet to train and validate ridge, lasso, and elastic net.
	1. Predict the phenotype value for both training and testing  folds using  λ  within 1 standard error of the minimum obtained by an inner 5-fold CV.  
	2. Record the appropriate performance metric using the optimal cutoff to optimize the metric.  
	3. Record the potentially significant SNPs from each method. The significant SNPs are those with coefficients higher than the mean of the absolute value of the coefficients.

	c. Filter SNPs by taking the union of SNPs from the ridge, LASSO, and elastic net.  
d. Create groups of SNPs using Hierarchical Clustering.  
e. Utilize filtered SNPs to train and validate Group Lasso and SGL using R functions grplasso and SGL.
	1. Predict the phenotype value for both training and testing folds using  λ  within 1 standard error of the minimum obtained by an inner 5-fold CV.  
	2. Record the appropriate performance metric using the optimal cutoff to optimize the metric.  
	3. Take the union of the potentially significant SNPs from both Group Lasso and SGL. The significant SNPs are those with coefficients higher than a cutoff (mean of the absolute value of the coefficients).

<b>Output:</b>  The significant SNPs (union of selected SNPs from Group LASSO and SGL) for each phenotype
