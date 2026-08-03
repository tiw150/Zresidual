# Acute myeloid leukemia survival data

Survival data for 1,043 adult acute myeloid leukemia patients in
northwest England, including residence coordinates, administrative
district, and subject-specific prognostic factors.

## Usage

``` r
data(LeukSurv)
```

## Format

A data frame with 1,043 rows and 9 variables:

- `time`:

  Survival time, in days.

- `cens`:

  Right-censoring status: 0 for censored and 1 for dead.

- `xcoord`:

  Residence coordinate on the x-axis.

- `ycoord`:

  Residence coordinate on the y-axis.

- `age`:

  Age in years.

- `sex`:

  Sex: 0 for female and 1 for male.

- `wbc`:

  White blood cell count at diagnosis, truncated at 500.

- `tpi`:

  Townsend deprivation score; higher values indicate less affluent
  areas.

- `district`:

  Administrative district of residence, among 24 districts.

## Source

Adapted from the `LeukSurv` data distributed in the spBayesSurv package.

## References

Henderson, R., Shimakura, S., and Gorst, D. (2002). Modeling spatial
variation in leukemia survival data. *Journal of the American
Statistical Association*, 97, 965-972.
[doi:10.1198/016214502388618753](https://doi.org/10.1198/016214502388618753)
.
