# Print methods for Z-residual objects

Compact print methods for `"zresid"` and `"cvzresid"` objects.

## Usage

``` r
# S3 method for class 'zresid'
print(x, ...)

# S3 method for class 'cvzresid'
print(x, ...)
```

## Arguments

- x:

  An object returned by
  [`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md)
  or
  [`Zresidual_CV`](https://tiw150.github.io/Zresidual/reference/Zresidual_CV.md).

- ...:

  Further arguments passed to
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The input object, invisibly.

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  set.seed(1)
  n <- 20
  x <- rnorm(n)
  t_event <- rexp(n, rate = exp(0.2 * x))
  t_cens  <- rexp(n, rate = 0.5)
  status  <- as.integer(t_event <= t_cens)
  time    <- pmin(t_event, t_cens)
  dat <- data.frame(time = time, status = status, x = x)
  fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)

  z <- Zresidual(fit, data=dat, nrep = 1, seed = 1)
  print(z)
}
#> Warning: NaNs produced
#> <zresid> n=20, reps=1
#>             rep1
#> [1,]  0.12371351
#> [2,] -0.67041677
#> [3,] -0.19062811
#> [4,] -0.80128050
#> [5,] -1.28481542
#> [6,] -0.05743129
```
