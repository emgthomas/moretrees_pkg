# moretrees

![](https://travis-ci.org/IQSS/moretrees.svg?branch=master) [![Github All Releases](https://img.shields.io/github/downloads/IQSS/moretrees/total.svg)]()

This package fits multi-outcome regression models to matched case-control and case-crossover data using Multi-Outcome Regression
with Tree-structured Shrinkage (MOReTreeS).

## Getting Started

MOReTreeS is a statistical method for analyzing the effect of an exposure on a large number of related outcomes.
The model performs three important functions simultaneously: (1) discovers groups of outcomes that have a similar relationship with the exposure;
(2) estimates group-specific exposure effects; (3) shares information between related outcomes about exposure effects and, where
relevant, covariate effects.
See the mortrees vignette for more information.

```{r}
vignette("moretrees")
```

If using moretrees in published research, please cite the following paper:

    Thomas EG, Trippa L, Parmigiani G, and Dominici F. Estimating the effects of fine particulate matter on 432 cardiovascular diseases using multi-outcome regression with tree-structured shrinkage. *Journal of the American Statistical Association*, 2020.

Currently, the package is only designed to fit MOReTreeS models to data that would normally be analyzed using conditional logistic regression, such as matched case-control or case-crossover data.
However, MOReTreeS equivalents of other regression models are in development.

### Installing

Moretrees can be installed by running the following code in your R session.
Use the option build_vignettes = TRUE to access the moretrees vignette.

```{r}
# install.packages("devtools")
library(devtools)
install_github("emgthomas/moretrees_pkg", build_vignettes = TRUE)
library(moretrees)
```

## Authors

Emma Thomas.
