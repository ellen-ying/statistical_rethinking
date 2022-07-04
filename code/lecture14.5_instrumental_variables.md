Lecture 14.5: Instrumental Variables
================
Yurun (Ellen) Ying
2022-07-04

## Instrumental variables

Correlated varying effects can also be implemented in constructing
instrumental variables for causal inferences.

Instrumental variables have to satisfy the following three criteria:

-   Independent of confounder U
-   Not independent of exposure X
-   Cannot influence outcome Y except through X

An example is th effect of education on wages. There are many
confounders, and we choose the season of birth as the instrumental
variable.

``` r
dag1 <- dagitty("dag{U[unobserved]; Q -> E -> W; E <- U-> W}")
coordinates(dag1) <- list(x = c(Q = 0.3, U = 0.5, E = 0.4, W = 0.6),
                          y = c(Q = 0.25, U = 0.25, E = 0.75, W = 0.75))
drawdag(dag1)
```

![](lecture14.5_instrumental_variables_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

A generative model of the data

``` r
set.seed(73)
N <- 500
U <- rnorm(N)
Q <- sample(1:4, size = N, replace = TRUE)
E <- rnorm(N, U + Q)
W <- rnorm(N, U) # education doesn't affect wages
dat_sim <- list(W = standardize(W),
                E = standardize(E),
                Q = standardize(Q))
```

The way to model this is just to write a statistical version of
generative model. To capture the effect of the unobserved confounders,
we will model E and W as from the same multivariate normal distribution:

![\begin{aligned}
\begin{pmatrix}
W_i\\\\
E_i
\end{pmatrix} &\sim \mathrm{MVNormal}(
\begin{pmatrix}
\mu\_{W,i}\\\\
\mu\_{E,i}
\end{pmatrix}, 
\mathbf{S}
) \\\\
\mu\_{W,i} &= \alpha_W + \beta\_{EW} E_i \\\\
\mu\_{E,i} &= \alpha_E + \beta\_{QE}Q_i
\end{aligned}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbegin%7Baligned%7D%0A%5Cbegin%7Bpmatrix%7D%0AW_i%5C%5C%0AE_i%0A%5Cend%7Bpmatrix%7D%20%26%5Csim%20%5Cmathrm%7BMVNormal%7D%28%0A%5Cbegin%7Bpmatrix%7D%0A%5Cmu_%7BW%2Ci%7D%5C%5C%0A%5Cmu_%7BE%2Ci%7D%0A%5Cend%7Bpmatrix%7D%2C%20%0A%5Cmathbf%7BS%7D%0A%29%20%5C%5C%0A%5Cmu_%7BW%2Ci%7D%20%26%3D%20%5Calpha_W%20%2B%20%5Cbeta_%7BEW%7D%20E_i%20%5C%5C%0A%5Cmu_%7BE%2Ci%7D%20%26%3D%20%5Calpha_E%20%2B%20%5Cbeta_%7BQE%7DQ_i%0A%5Cend%7Baligned%7D "\begin{aligned}
\begin{pmatrix}
W_i\\
E_i
\end{pmatrix} &\sim \mathrm{MVNormal}(
\begin{pmatrix}
\mu_{W,i}\\
\mu_{E,i}
\end{pmatrix}, 
\mathbf{S}
) \\
\mu_{W,i} &= \alpha_W + \beta_{EW} E_i \\
\mu_{E,i} &= \alpha_E + \beta_{QE}Q_i
\end{aligned}")

``` r
m <- ulam(
  alist(
    c(W,E) ~ multi_normal(c(muW,muE), Rho, Sigma),
    muW <- aW + bEW*E,
    muE <- aE + bQE*Q,
    c(aW,aE) ~ normal(0, 0.2),
    c(bEW,bQE) ~ normal(0, 0.5),
    Rho ~ lkj_corr(2),
    Sigma ~ exponential(1)
  ), data = dat_sim, chains = 4, cores = 4
)
```

``` r
round(precis(m, depth = 3), 2)
```

    ##           mean   sd  5.5% 94.5%   n_eff Rhat4
    ## aE        0.00 0.03 -0.06  0.05 1287.37     1
    ## aW        0.00 0.05 -0.08  0.07 1441.71     1
    ## bQE       0.59 0.04  0.53  0.65 1312.64     1
    ## bEW      -0.04 0.08 -0.17  0.07 1018.13     1
    ## Rho[1,1]  1.00 0.00  1.00  1.00     NaN   NaN
    ## Rho[1,2]  0.54 0.05  0.45  0.62  809.52     1
    ## Rho[2,1]  0.54 0.05  0.45  0.62  809.52     1
    ## Rho[2,2]  1.00 0.00  1.00  1.00     NaN   NaN
    ## Sigma[1]  1.02 0.05  0.95  1.10 1058.94     1
    ## Sigma[2]  0.81 0.02  0.77  0.85 1749.66     1

The influence of E on W is around zero, which reflects the true causal
influence. The correlation between E and W is reliably positive,
reflecting the common cause U.
