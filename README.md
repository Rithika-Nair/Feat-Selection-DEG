# Exploring the Impact of Knowledge-Driven Feature Selection on Human Differential Expression Analysis

## Abstract
Differential expression (DE) analysis is a crucial computational method for identifying genes that exhibit varying expression levels between different conditions, aiding in the understanding of biological processes and disease mechanisms. However, the high dimensionality and complexity of gene expression data pose challenges such as overfitting in machine learning models. This study investigates the impact of incorporating biological knowledge into feature selection methods, specifically focusing on a gene’s prior probability of differential expression, as well as mRNA and protein half-lives. Using data from The Cancer Genome Atlas across four cancer types, we implemented various feature selection strategies and evaluated their effects on model performance through cross-validated elastic net modeling. Our findings reveal that while the DE prior did not enhance predictive performance, the stability measures of mRNA and protein half-lives showed cancer-type-specific improvements in classification accuracy. These results underscore the importance of integrating biological knowledge into feature selection processes, suggesting that certain features can significantly influence model quality in gene expression analysis.



