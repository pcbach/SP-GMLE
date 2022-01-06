# covariance-constrained Gaussian Maximum Likelihood Estimation via semidefinite approximations (SP-GMLE)

**James Saunderson** James.Saunderson@monash.edu

**Chi Pham** Chi.Pham@monash.edu

This is a R parser for CVXR to solve the Gaussian MLE problem with added constraints. 

### Update

 - **Version 1.0** Initial release

### Table of content

- Uses
- Setup
- Basic uses
- Constraint
- Debug

## Introduction

This parser is use to solve the [MLE](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) problem for [p-dimensional Gaussian models](https://en.wikipedia.org/wiki/Multivariate_normal_distribution) with convex constraints on the covariance matrix. Which is a generalization of the classical linear covariance models. 

The problem was proven to be convex if there are enough data by [Zwiernik, Uhler, and Richards](https://arxiv.org/pdf/1408.5604.pdf). This parser produce a epigraph of the negative-log likelihood function by using [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature), a weighted-sum of positive semi-definite epigraphs.

## Setup
To use this parser, you need an installation of [R](https://www.rstudio.com/) as well as [CVXR](https://cvxr.rbind.io/). Other than that, [MOSEK](https://www.mosek.com/) and [Rmosek](https://cran.r-project.org/web/packages/Rmosek/index.html) is highly recommended if you want to add any constraints.

```
rm(list = ls())
suppressMessages(require("Rmosek"))
suppressMessages(library("CVXR"))
suppressMessages(source("SPGL.R"))
#set.seed(124)
msk = MOSEK()
import_solver(msk)
```

## Basic uses
The parser take in 5 parameters:
- **d**: problem dimension
- **Sigma**: dxd CVXR variable
- **Tau**: CVXR scalar variable
- **Sn**: [sample covariance](https://en.wikipedia.org/wiki/Sample_mean_and_covariance#Unbiasedness) of the data
- **n**: Gaussian quadrature degree

```

#randomly generated data
data = rnorm(4*20) #20 datapoints
data = matrix(data,nrow = 20)
Sn = cov(data,data)*(dim(data)[1]-1)/(dim(data)[1])
#rescaled from unbiased estimator

n = 5
d = dim(data)[2]

Sigma <- Variable(d,d)
Tau <- Variable(1)

#parser generate epigraph constraints
rep <- represent(Sigma,Tau,n,d,Sn)
constraints <- rep$constraints

#normal CVXR syntax
objective <- Minimize(Tau)
problem <- Problem(objective, constraints)
result <- solve(problem, solver = "MOSEK")

print(result$getValue(Sigma))
print(result$getValue(Tau))
```

Any CVXR variable define outside the parser (Sigma and Tau in this case) can be access by using ```result$getValue(variableName)```

## Constraint
To add constraint on the covariance matrix, the CVXR variable Sigma can be conditioned directly. Below is an example of a correlation matrix constraint (diagonal entry equal 1)
```
constraints <- c(constraints,list(
  diag(Sigma) == rep(1,d)
))
```

There are limitations on which way the constraints can be add, most of the time you will need to use the [inbuilt CVXR functions](https://cvxr.rbind.io/cvxr_functions/). Below is an example of a more complicated constraint, a brownian motion model.
```
#Tree description matrix
E <- matrix(
  c(1, 0, 0, 0, 
    0, 1, 0, 0, 
    0, 0, 1, 0,
    0, 0, 0, 1,
    1, 1, 1, 1
    ),  
  ncol = d,        
  byrow = TRUE ) 
  
#Extra variables
r = dim(E)[1]  # Number of nodes
v <- Variable(r, nonneg=TRUE) #Branch length

#Brownian motion tree constraint
for (i in seq(from = 1, to = d,by = 1)){
    for (j in seq(from = 1, to = d,by = 1)){
        ancestor = rep(0,r)
        #Boolean, ancestor[k] = 1 if k is an common ancestor of i and j   
        for (k in seq(from = 1, to = r,by = 1)){
            EtE = E[k,] %*% t(E[k,])
            if(EtE[i,j] == 1){
                ancestor[k] = 1
            }
        }
        #dot product a.v
        temp = sum_entries(ancestor*v)
        constraints <- c(constraints,list(
            X[i,j] == temp
        ))
    }
}


print(result$getValue(v))
```

## Debug
In cases where you want to step deep into the parser itself, there are options to access the underlying variables, the ids of these underlying variable is stored in ```rep$variableid```. For example to access variable "x", ```result[[ids[["x"]]]]```.

