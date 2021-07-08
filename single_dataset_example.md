Example of running the functions with a single simulated dataset
================
M Schnitzer
08/07/2021

## Generate a single data set

We generate a single simulated data set using the function datagen in
the file `datagen_func.R`.

``` r
source("datagen_func.R")
seed<-20210708
```

We must set three parameters:

-   `beta`, the strength of the effect of *D*(2) and *D*(3) on *Y* (note
    that beta is on the odds ratio scale in the paper but here it is on
    the log-odds scale)
-   `theta1`, the strength of the effect of *A*(2) on *Y*, log-odds
    scale (dashed arrow in Figure 3)
-   `theta2`, the strength of the effect of *A*(4) on *Y*, log-odds
    scale (dotted arrow in Figure 3)

These parameters are set as follows in the paper’s simulation study:

-   scenario 1: `beta=theta1=theta2=0`,
-   scenario 2: `beta=log(2); theta1=theta2=0`,
-   scenario 3: `beta=log(2); theta1=0.4; theta2=0` – only scenario
    where there is an effect of *A*(2) on *Y*,
-   scenario 4: `beta=log(2); theta1=0; theta2=0.4` – only scenario
    where there is an effect of *A*(4) on *Y*,

``` r
beta=log(2); theta1=0; theta2=0.4

DAT<-datagen(ssize=1000,seed=seed,beta=beta,theta1=theta1,theta2=theta2)
```

`DAT` is a data frame that has four exposure nodes. The variable
ordering in the data frame is assumed to be
(*W*(1), *A*(1), *W*(2), *A*(2), *D*(2), *W*(3), *A*(3), *D*(3), *W*(4), *A*(4), *Y*).

Define the desired initiation time *k* in {2, 3, 4, 5} where 5 means
never initiate. Here we set *k* = 3 (initiate at time 3) which gives the
exposure pattern (0, 0, 1, 1).

``` r
k<-3
```

## Apply the four estimators

First we load the estimator functions from the file `estimators_func.R`.

``` r
source("estimators_func.R")
```

Then we can run all target trials methods, which output the estimation
of *E*(*Y*<sup>*k* = 3</sup>).

``` r
ipwfunc(dat=DAT,k=k)
```

    ## (Intercept) 
    ##    0.234569

``` r
gcompfunc(dat=DAT,k=k)
```

    ## [1] 0.2470772

``` r
tmlefunc(dat=DAT,k=k)
```

    ## [1] 0.2356928

I also included the \`\`standard" IPW method that ignores delivery time
but otherwise adjusts for time-dependent confounding.

``` r
ipw_standfunc(dat=DAT,k=k)
```

    ## (Intercept) 
    ##   0.1614631

## Approximating the true value of *E*(*Y*<sup>*k* = 3</sup>)

What is the true value of *E*(*Y*<sup>*k* = 3</sup>) under this data
generation? The `datagen` function also includes an option to generate
counterfactual data under a given treatment initiation time, *k*.
Without setting the seed, it will generate a random seed.

``` r
DAT_cf<-datagen(ssize=100000,beta=beta,theta1=theta1,theta2=theta2,counterfactual=T,k=k)
```

    ## Generating a counterfactual data set with treatment pattern 0 0 1 1

``` r
mean(DAT_cf$Y)
```

    ## [1] 0.24766

You can run this with a greater sample size and/or multiple times to get
good approximations of the true value.
