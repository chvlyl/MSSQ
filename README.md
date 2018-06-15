# MSSQ
Multi-Sample Poisson Model for Species Abundance Quantification

## Introduction
The human microbiome, which includes the collective microbes residing in or on the human body, has a profound influence on the human health. DNA sequencing technology has made the large-scale human microbiome studies possible by using shotgun metagenomic sequencing. One important aspect of data analysis of such metagenomic data is to quantify the bacterial abundances based on the metagenomic sequencing data. Existing methods almost always quantify such abundances one sample at a time, which ignore certain systematic differences in read coverage along the genomes due to GC contents, copy number variation and the bacterial origin of replication. 

In order to account for such differences in read counts, we propose a multi-sample Poisson model to quantify microbial abundances based on read counts that are assigned to species-specific taxonomic markers. Our model takes into account the marker-specific effects when normalizing the sequencing count data in order to obtain more accurate quantification of the species abundances. Compared to currently available methods on simulated data and real data sets, our method has demonstrated an improved accuracy in bacterial abun- dance quantification, which leads to more biologically interesting results from downstream data analysis.


## Installation
You can install our MSSQ package from Github
```r
install.packages("devtools")
devtools::install_github("chvlyl/MSSQ")
library(MSSQ)
```

## Basic Usage
```r
est <- poisson_estimate_parameter_LS(X=X,tc=tc,l=l,N=N,S=S,
                                  species.names=species.names,
                                  estimate.phi=TRUE)
```

## Details
Consider a metagenomic study with N samples. After the sequencing reads are aligned to
sets of clade-specific marker genes, the data can be summarized as a large table of counts, where X_ijk is the count data of sequencing reads for sample i (i = 1, 2, …, n), species j (j = 1, 2, …, p) and marker k (k = 1, 2, …, m_j). We model the count data for all species and all samples together and assume that the count X_ijk is generated from the following Poisson model,

<p align="center">
  <img src="images/eqn1.png" width="350"/>
</p>

<p align="center">
  <img src="images/eqn2.png" width="350"/>
</p>




## Citation
Eric Z. Chen, Frederic D. Bushman and Hongzhe Li (2015). A Model-Based Approach For Species Abundance Quantification Based On Shotgun Metagenomic Data. Submitted.
