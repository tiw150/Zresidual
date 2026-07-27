# Inverse Gaussian density

Density function for the inverse Gaussian distribution using the same
parameterization as brms' inverse.gaussian family: mean parameter `mu`
and shape parameter `shape`.

## Usage

``` r
dinvgauss(x, mu, shape, log = FALSE)
```

## Arguments

- x:

  Numeric vector of positive values.

- mu:

  Positive mean parameter.

- shape:

  Positive shape parameter.

- log:

  Logical; if `TRUE`, return log-density.

## Value

A numeric vector of density values or log-density values.
