# SGM
An R package for Estimating Spectral Graphical Models

*** 
The R package `SGM` provides a set of functions that learn a latent graph from the data.  The data are viewed as stationary signals from the inferred latent graph.  `SGM` is motivated by the need to understand complex dependence structures, and how they change over time.  Currently, the package offers code to fit a GAR(1) model (as proposed by Peng et al.) through a 3-step procedure.

*** 

## Dependencies 
Please make sure that the following packages are installed before using the R package `SGM`. 

```
install.packages(c("doParallel", "foreach", "gmp"))
```

## Installation
The R Package `SGM` can be installed from source files, in this GitHub repository.  The `devtools` package is required for installation.  
```
library(devtools)
install_github(repo="jed-harwood/SGM")
```

## Main Functions
`GAR1_ADMM_fit`: learn a GAR(1) SGM for a given set of tuning parameters

`model_selection`: conduct model selection via the eBIC criterion or log-likelihood

## Contact
Please report any bugs to `jedharwood@ucdavis.edu`
