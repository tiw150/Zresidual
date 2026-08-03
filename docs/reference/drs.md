# Diabetic Retinopathy Study data

Paired-eye survival data from 197 high-risk diabetic patients. One
randomly selected eye received laser treatment and the other eye was
untreated; blindness was the event of interest.

## Usage

``` r
data(drs)
```

## Format

A data frame with 394 rows and 8 variables, with two rows for each
patient:

- `subject_id`:

  Unique patient identifier.

- `eye`:

  Eye indicator: 1 for right and 2 for left.

- `time`:

  Observed follow-up time.

- `status`:

  Outcome indicator: 0 for censored and 1 for blindness.

- `treated`:

  Treatment indicator: 0 for untreated and 1 for treated.

- `age_at_onset`:

  Age at onset of diabetes, in years.

- `laser_type`:

  Laser type: 1 for xenon and 2 for argon.

- `diabetes_type`:

  Diabetes type: 1 for juvenile and 2 for adult.

## Source

Adapted from the `drs` data distributed in the frailtySurv package.

## References

Huster, W. J., Brookmeyer, R., and Self, S. G. (1989). Modelling paired
survival data with covariates. *Biometrics*, 45, 145-156.
[doi:10.2307/2532041](https://doi.org/10.2307/2532041) .
