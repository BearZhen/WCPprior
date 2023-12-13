# WCPprior - Wasserstein complexity penalization priors

[![R-CMD-check](https://github.com/vpnsctl/WCPprior/actions/workflows/r-check-devel.yml/badge.svg)](https://github.com/vpnsctl/WCPprior/actions/workflows/r-check-devel.yml)

## Overview 

`WCPprior` is an R package used to obtain numerical approximation of Wasserstein complexity penalization (WCP) priors. The package also contains several analytical WCP priors as well as implementations for STAN and INLA.

Please visit the [package homepage][ref2] for detailed tutorials on all different aspects of the package.

# Installation instructions #

The latest development version can be installed by using the command
```r
remotes::install_github("vpnsctl/WCPprior", ref = "devel")
```

The latest stable version can be installed by using the command
```r
remotes::install_github("vpnsctl/WCPprior", ref = "stable")
```

# References #
D. Bolin, A. Simas, Z. Xiong (2023+) [Wasserstein complexity penalization priors: a new class of penalizing complexity priors][ref]. arXiv preprint.



[ref]: https://arxiv.org/abs/2312.04481  "Wasserstein complexity penalization priors: a new class of penalizing complexity priors"

[ref2]: https://vpnsctl.github.io/WCPprior/ "WCPprior homepage"