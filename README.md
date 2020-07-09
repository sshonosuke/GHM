# GHM: Grouped Heterogeneous Mixture Modeling for Clustered Data
This package implements mixture models for clustered data, as proposed by the following paper.

Sugasawa, S. (2020). Grouped Heterogeneous Mixture Modeling for Clustered Data. *Journal of the American Statistical Association*, to appear.  https://arxiv.org/abs/1804.00888

Functions are implemented in *GHM-Gauss-function.R* and *GHM-Po-function.R* for continous and count responses, respectively.


# Continuous response 
Install function and demo data.
```{r}
source("GHM-Gauss-function.R")
load("simdata-Gauss.RData")
```

Two functions are included in *GHM-Gauss-function.R*.
`GHM.Gauss` fits the proposed method with specified numbers of groups `G` and latent distributions `L`.
`aGHM.Gauss` fits the adaptive version of the proposed method with `G` and `L` selected via an information criterion.

Apply the adaptive version of the proposed method. 

Input of `aGHM.Gauss`
- `Y`: vector of continous response variables  
- `X`: (N,p)-matrix of covariates without intercept term (N: total sample size; p: the number of covariates)
- `ID`: vector of cluster ID 
- `rG`: range of candidates for the number of groups `G`
- `rL`: range of candidates for the number of latent components `L`
- `print`: If `T`, progress will be output

Output of `aGHM.Gauss`
- `fit.best`: result of the best model
- `fit.all`: list of results of all the candidate models 
- `BIC`: matrix of values of the information criterion 

`fit.best` is a list object including the following information
- `ML`: maximum likelihood value
- `Beta`: matrix of estimated regression coefficients for each component  
- `Sig`: vector of estimated standard deviations 
- `Pi`: matrix of estimated mixing proportions for each group
- `g`: vector of estimated grouping parameters
- `BIC`: value of the information criterion 

```{r}
fit <- aGHM.Gauss(Y, X, ID, rG=c(2, 10), rL=c(2, 4), print=T)
fit$BIC     # information criterion for all the candidate models

bfit <- fit$fit.best
bfit$Beta   # regression coefficients
bfit$g      # grouping parameters
bfit$Pi     # group-wise mixing proportions
```


# Count response 
Install function and demo data.
```{r}
source("GHM-Po-function.R")
load("simdata-Po.RData")
```

Two functions `GHM.Gauss` (fixed `G` and `L`) and `aGHM.Gauss` (adaptive version) are included in *GHM-Po-function.R*.
Inputs and outputs of these functions are the same as the continuous case.　 

```{r}
fit <- aGHM.Po(Y, X, ID, rG=c(2, 10), rL=c(2, 4), print=T)
fit$BIC     # information criterion for all the candidate models

bfit <- fit$fit.best
bfit$Beta   # regression coefficients
bfit$g      # grouping parameters
bfit$Pi     # group-wise mixing proportions
```
 
