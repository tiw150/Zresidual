# A function to calculate Bartlett of Zresidual

A function to calculate Bartlett of Zresidual

## Usage

``` r
# S3 method for class 'zresid'
bartlett.test(x, X = c("lp", "covariate"), k.bl = 10, ...)
```

## Arguments

- x:

  A Z-residual object (with class 'zresid').

- X:

  Linear predictor or covariate. Must be (1) a length-n vector, (2)
  'lp'/'covariate', or (3) a covariate name in attr(x, 'covariates').

- k.bl:

  Number of bins if applicable for continuous covariates. Default is 10.

- ...:

  Further arguments passed to or from other methods.
