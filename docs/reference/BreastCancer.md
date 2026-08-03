# German Breast Cancer Study Group 2 data

A numeric recoding of recurrence-free survival data for 686 women with
node-positive breast cancer in the German Breast Cancer Study Group 2
trial.

## Usage

``` r
data(BreastCancer)
```

## Format

A data frame with 686 rows and 12 variables:

- `X`:

  Row index in the supplied data file.

- `id`:

  Patient identifier.

- `time`:

  Recurrence-free survival time, in days.

- `status`:

  Event indicator: 0 for censored and 1 for recurrence or death.

- `treat`:

  Tamoxifen indicator: 0 for no tamoxifen and 1 for tamoxifen.

- `age`:

  Age in years.

- `men`:

  Menopausal status: 1 for premenopausal and 2 for postmenopausal.

- `size`:

  Tumor size, in millimeters.

- `grade`:

  Tumor grade, coded 1 to 3.

- `nodes`:

  Number of positive lymph nodes.

- `prog`:

  Progesterone receptor level, in fmol.

- `oest`:

  Estrogen receptor level, in fmol.

## Source

Adapted from the `GBSG2` data distributed in the TH.data package.

## References

Schumacher, M., et al. (1994). Randomized 2 x 2 trial evaluating
hormonal treatment and the duration of chemotherapy in node-positive
breast cancer patients. *Journal of Clinical Oncology*, 12, 2086-2093.
[doi:10.1200/JCO.1994.12.10.2086](https://doi.org/10.1200/JCO.1994.12.10.2086)
.
