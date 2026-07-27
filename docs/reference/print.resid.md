# Print a residual vector

Simple print method for residual-like vectors used in this package.

## Usage

``` r
# S3 method for class 'resid'
print(x, ...)
```

## Arguments

- x:

  A residual object or numeric vector.

- ...:

  Further arguments passed to
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The input object, invisibly.

## Examples

``` r
x <- structure(rnorm(5), class = "resid")
print(x)
#> [1]  0.9969869 -0.2757780  1.2560188  0.6466744  1.2993123
```
