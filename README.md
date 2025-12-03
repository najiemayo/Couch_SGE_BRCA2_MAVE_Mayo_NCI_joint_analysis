# Couch_SGE_BRCA2_MAVE_Mayo_NCI_joint_analysis

This repository provides the VarCall model tool used in the analysis of high throughput CRISPR based saturated genome editing (SGE) data, one type of multiplexed assays of variant effects (MAVEs). It includes two data sets joint analysis. The VarCall model mainly uses R and the R package rjags to perform the Bayesian two Gaussian components modeling.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Test and Run](#Test-and-run)
- [License](#license)

# Overview
``VarCall`` aims to provide a comprehensive statistical model for SGE data analysis, including the modeling of batch effects, and the prediction of variant effects. The package utilizes a Bayesian hierarchical two Gaussian components modeling. The package should be able to run on all major platforms (e.g. BSD, GNU/Linux, OS X, Windows) as long as R and related R Packages are installed.

# System Requirements
## Hardware requirements
`VarCall` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for *Linux*, but in principle should be able to run on all other platforms such as OS X and Windows. The package has been tested on the following systems:
+ Linux: Red Hat Enterprise Linux 8.8
+ Windows: Windows 10 Enterprise

### Package Dependencies
`VarCall` mainly depends on JAGS and the following R packages.

```
knitr, rjags, R2WinBUGS, R2jags, mgcv
```

# Installation Guide:
### Download or clone from github
```
git clone https://github.com/najiemayo/Couch_SGE_BRCA2_MAVE_Mayo_NCI_joint_analysis/VarCall/ 
```

### Install packages
First install JAGS in the system following the link here: https://mcmc-jags.sourceforge.io/
Then within R, type the following:
```
install.packages(c("knitr", "rjags", "R2WinBUGS", "R2jags", "mgcv"))
```
This should be done within one miniute.

# Test and Run:


- To run :
  - Download the full data files `combined.raw.tsv`, `variant_type_for_train.csv`, `shyamBRCA2allCounts.txt` and put them in the VarCall folder
  - In terminal, run R script. "Rscript runRtex.R"
  - Check if your working directory corresponds to where the script is looking for. 

- Specified prior:
  - Currently the prior is set to a mean value of 0.2 using a beta distribution Beta(2, 8). The change the prior, in `BRCA2CombinedMave24.ldaER.Run1.Rtex` file, modifiy the line `beta.a<-2.0` and `beta.b<-8.0`. 

- Expected output and running time.

File `MAVEpostProbs.csv` is the main output file for predicted probabilities of being pathogenic based on the training labels from file `variant_type_for_train.csv`. The following shows the main columns of output:
  - PrDel: the probability of being pathogenic
  - lPostOdds: the log posterior odds
  - logBF: the log Bayes factor
  - eta: the estimated effect size
  - eta.ll, eta.ul: 95% lower bound and upper bound of eta

Visualization of intermediate results:
  - lPOrun1vs2.pdf
  - PrDelHistogram.pdf
  - PrDelRun1vs2.pdf 
  - ResidualQQ.pdf
  - BatchLocationVsScale.pdf
  - BatchMeanBPlot.pdf
  - SmoothFitsFnl.pdf
  - SmoothFitsByRepFnl.pdf
  - BatchMeanQQ.pdf
  - VariantQQ.pdf           
  - BatchScaleBPlot.pdf 
  - BatchScaleQQ.pdf 
  - EtaHistogram.pdf
  - lPObarplots.pdf

Intermediate files:
  - mayoPreprocessing.tex
  - nciPreprocessing.tex
  - BRCA2CombinedMave24.ldaER.Run1.tex         
  - BRCA2CombinedMave24.ldaER.Run2.tex  
  - mayoData.RData
  - nciPreprocessing.Rtex
  - nciDB.RData
  - run1.RData

Two intermediate folders will be generated:
  - cache                   
  - figs

change setting of running `mcmc.pars` in `BRCA2CombinedMave24.ldaER.Run1.Rtex` file. 

%% begin.rcode\
%  ## T degrees of freedom for measurement error model\
%  t.df<-5\
%  ##\
%  ## Prior Probability Pathogenic, beta distribution parameters:\
%  beta.a<-2.0\
%  beta.b<-8.0  ##mean=0.20, ESS = 10\
%  ## Set List of MCMC Control Parameters:\
%  mcmc.pars<-list(iter=10000, ## short run\
%                  burn=5000,\
%                  thin=10)\
%  ##mcmc.pars<-list(iter=150000, ## long run\
%  ##               burn=50000,\
%  ##               thin=10)\
%  ##mcmc.pars<-list(iter=550000, ## longer run\
%  ##                 burn=50000,\
%  ##                 thin=25)\
%% end.rcode\

Three settings with parameters are prepared in the file, comment with `##` or uncomment by removing `##` for `mcmc.pars`.

The running time for the short run takes about one hour. The running time for the longer run takes ~4days. Due to the stochastic nature of the model, there might be some minor differences in the results, but it should converge when the iteration is long enough.

# License

This project is covered under the **MIT License**.


