# m6Aexpress-enet
## Introduction
As the most abundant mRNA modification, m6A controls and influences many aspects of mRNA metabolism including the mRNA stability and degradation. However, the role of specific m6A sites in regulating gene expression still remains unclear. In additional, the multicollinearity problem caused by the correlation of methylation level of multiple m6A sites in each gene could influence the prediction performance. To address the above challenges, we propose an enet-regularization negative binomial regression model (called m6Aexpress-enet) to predict which m6A site could potentially regulate its gene expression. 
# The step-by-step implementation of m6Aexpress-enet
## Installed the Required R packages
Before using m6Aexpress-enet, we should firstly install some required R packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsamtools","GenomicAlignments","GenomicRanges",
                       "GenomicFeatures","rtracklayer","exomePeak2"))
install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/DESeq_1.34.1.tar.gz", repos = NULL, type="source")
install.packages("randomForest")
```
## Dowloand SRAMP for predicting m6A sites in single base resolution
1. Go to the webserver: https://www.cuilab.cn/sramp/ and download SRAMP tool

2. Require the Perl (version>5.8)



