# Kidney catheter infection data

Recurrence times to infection at the catheter insertion point for 38
kidney patients using portable dialysis equipment. Each patient has
exactly two observations.

## Usage

``` r
data(kidney)
```

## Format

A data frame with 76 rows and 7 variables:

- `id`:

  Patient identifier.

- `time`:

  Observed recurrence or censoring time.

- `status`:

  Event indicator: 0 for censored and 1 for infection.

- `age`:

  Age in years.

- `sex`:

  Sex: 1 for male and 2 for female.

- `disease`:

  Disease type, a factor with levels `Other`, `GN`, `AN`, and `PKD`.

- `frail`:

  Frailty estimate reported in the original study.

## Source

Adapted from the `kidney` data distributed in the survival package.

## References

McGilchrist, C. A. and Aisbett, C. W. (1991). Regression with frailty in
survival analysis. *Biometrics*, 47, 461-466.
[doi:10.2307/2532138](https://doi.org/10.2307/2532138) .
