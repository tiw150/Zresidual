# Z-residual diagnostic tool for assessing covariate functional form in shared frailty models

## 1 Installing Zresidual and Other packages

### 1.1 Installing Z-residua from the source

Code

``` r
devtools::install()
```

    ── R CMD build ─────────────────────────────────────────────────────────────────
    * checking for file ‘/Users/longhai/Github/Zresidual/DESCRIPTION’ ... OK
    * preparing ‘Zresidual’:
    * checking DESCRIPTION meta-information ... OK
    * checking for LF line-endings in source and make files and shell scripts
    * checking for empty or unneeded directories
    Removed empty directory ‘Zresidual/data’
    Removed empty directory ‘Zresidual/vignettes/_freeze’
    Omitted ‘LazyData’ from DESCRIPTION
    * building ‘Zresidual_0.1.0.tar.gz’

    Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
      /var/folders/mh/fxdzznn94tg5z23jqm1n3tgc0000gq/T//RtmpVgkmnf/Zresidual_0.1.0.tar.gz \
      --install-tests
    * installing to library ‘/Users/longhai/Library/R/arm64/4.5/library’
    * installing *source* package ‘Zresidual’ ...
    ** this is package ‘Zresidual’ version ‘0.1.0’
    ** using staged installation
    ** R
    ** inst
    ** byte-compile and prepare package for lazy loading
    ** help
    *** installing help indices
    *** copying figures
    ** building package indices
    ** installing vignettes
    ** testing if installed package can be loaded from temporary location
    ** testing if installed package can be loaded from final location
    ** testing if installed package keeps a record of temporary installation path
    * DONE (Zresidual)

Code

``` r
library(Zresidual)
```

Code

``` r
if (!requireNamespace("Zresidual", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("tiw150/Zresidual",
                          upgrade = "never",
                          dependencies = TRUE)
}
library(Zresidual)
```

### 1.2 Intalling and Loading R Packages used in this Demo

Code

``` r
pkgs <- c(
  "survival","EnvStats","foreach","statip","VGAM","plotrix","actuar",
  "stringr","Rlab","dplyr","rlang","tidyr",
  "matrixStats","timeDate","katex","gt","loo"
)

missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, dependencies = TRUE)
}

invisible(lapply(pkgs, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

nc <- parallel::detectCores(logical = FALSE)
if (!is.na(nc) && nc > 1) options(mc.cores = nc - 1)
```

## 2 Introduction

This vignette explains how to use the Zresidual package to calculate
Z-residuals based on the output of the coxph function from the survival
package in R. It also serves as a demonstration of how to use
Z-residuals to assess the overall goodness of fit (GOF) and identify
specific model misspecifications in semi-parametric shared frailty
models. To fully understand the detailed definitions and the example
data analysis results, please refer to the original paper titled
“Z-residual diagnostics for detecting misspecification of the functional
form of covariates for shared frailty models.

## 3 Definition of Z-residual

We use Z-residual to diagnose shared frailty models in a Cox
proportional hazard setting with a baseline function unspecified.
Suppose there are g groups of individuals, with each group containing
n_i individuals, indexed as i = 1, 2, , g in the case of clustered
failure survival data. Let y\_{ij} be a possibly right-censored
observation for the jth individual from the ith group, and \delta\_{ij}
be the indicator for being uncensored. The normalized randomized
survival probabilities (RSPs) for y\_{ij} in the shared frailty model is
defined as: \begin{equation} S\_{ij}^{R}(y\_{ij}, \delta\_{ij}, U\_{ij})
= \left\\ \begin{array}{rl} S\_{ij}(y\_{ij}), & \text{if \$y\_{ij}\$ is
uncensored, i.e., \$\delta\_{ij}=1\$,}\\ U\_{ij}\\S\_{ij}(y\_{ij}), &
\text{if \$y\_{ij}\$ is censored, i.e., \$\delta\_{ij}=0\$,} \end{array}
\right. \label{rsp} \end{equation} where U\_{ij} is a uniform random
number on (0, 1), and S\_{ij}(\cdot) is the postulated survival function
for t\_{ij} given x\_{ij}. S\_{ij}^{R}(y\_{ij}, \delta\_{ij}, U\_{ij})
is a random number between 0 and S\_{ij}(y\_{ij}) when y\_{ij} is
censored. It is proved that the RSPs are uniformly distributed on (0,1)
given x\_{i} under the true model . Therefore, the RSPs can be
transformed into residuals with any desired distribution. We prefer to
transform them with the normal quantile: \begin{equation}
r\_{ij}^{Z}(y\_{ij}, \delta\_{ij}, U\_{ij})=-\Phi^{-1}
(S\_{ij}^R(y\_{ij}, \delta\_{ij}, U\_{ij})),\label{zresid}
\end{equation} which is normally distributed under the true model, so
that we can conduct model diagnostics with Z-residuals for censored data
in the same way as conducting model diagnostics for a normal regression
model. There are a few advantages of transforming RSPs into Z-residuals.
First, the diagnostics methods for checking normal regression are rich
in the literature. Second, transforming RSPs into normal deviates
facilitates the identification of extremely small and large RSPs. The
frequency of such small RSPs may be too small to be highlighted by the
plots of RSPs. However, the presence of such extreme SPs, even very few,
is indicative of model misspecification. Normal transformation can
highlight such extreme RSPs.

## 4 Examples for Illustration and Demonstration

### 4.1 Load the real Dataset

This example provides a fundamental illustration of using the
Z-residuals for diagnosing both the overall goodness of fit (GOF) and
the functional form of covariates in a real application for modelling
the survival times of acute myeloid leukemia patients.

The dataset employed in our analysis contains 411 patients who were
recorded at the M. D. Anderson Cancer Center between 1980 and 1996.
These patients were under the age of 60 and hailed from 24
administrative districts. The data collected information on the survival
time for acute myeloid leukemia and prognostic factors, including age,
sex, white blood cell count (wbc) at diagnosis, and the townsend score
(tpi) for which higher values indicate less affluent areas. The
censoring rate is 29.2%. The response variable of interest is the
survival time in days, which is the time from entry to the study or
death. The preliminary study showed that the wbc is highly right-skewed.
Logarithm transformation is often used to reduce the impact of extremely
large values of the covariate on the response variable, such as the wbc
variable in this application. However, a logarithm transformation may
mask the impact of extremely large values of the covariate on the
outcome variable.

Code

``` r
data_path <- system.file("extdata", "LeukSurv.rda", package = "Zresidual")
load(data_path)

LeukSurv <- transform(LeukSurv,
  district = as.factor(district),
  sex      = as.factor(sex),
  logwbc   = log(wbc + 0.001)
)

LeukSurv <- LeukSurv[LeukSurv$age < 60, ]
```

### 4.2 Fitting Models

We fitted two shared frailty models, one with covariates wbc, age, sex
and tpi, which is labelled as the wbc model, and the other with log(wbc)
replacing wbc, which is labelled as the lwbc model.

Code

``` r
fit_LeukSurv_wbc <- coxph(Surv(time, cens) ~ age  +sex+ wbc +tpi  +
          frailty(district, distribution="gamma"), data= LeukSurv)
fit_LeukSurv_logwbc  <- coxph(Surv(time, cens) ~ age +sex + logwbc + tpi +
          frailty(district, distribution="gamma"), data= LeukSurv)
```

### 4.3 Computing Z-Residuals

Once the model is fitted, we can calculate Z-residuals for two models.

Code

``` r
Zresid.LeukSurv.wbc<-Zresidual(object = fit_LeukSurv_wbc,nrep=1000)
Zresid.LeukSurv.logwbc<-Zresidual(object = fit_LeukSurv_logwbc,nrep=1000)
```

### 4.4 Inspecting the Normality of Z-Residuals for Checking Overall GOF

Diagnosing the overall goodness-of-fit (GOF) using Z-residuals as
follows:

A QQ plot based on Z-residuals can be used to graphically assess the
model’s overall GOF, and Shapiro-Wilk (SW) or Shapiro-Francia (SF)
normality tests applied to Z-residuals can be used to numerically test
the overall GOF of the model. We can see that the QQ plots of
Z-residuals of these two models align well with the 45 ^\circ diagonal
line. The Z-SW tests also give large p-values for two models, where Z-SW
is the test method that the normality of Z-residuals is tested with the
SW test.

Code

``` r
for (i in 1:10) {
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
   qqnorm.zresid(Zresid.LeukSurv.wbc,irep=i)
   qqnorm.zresid(Zresid.LeukSurv.logwbc, irep=i)
}
```

![](zresidual_demo_files/figure-html/qqplot-zresid-anim-.gif)

Figure 1: QQ plots of Z-residuals for the wbc (left panels) and lwbc
(right panels) models fitted to the survival data of acute myeloid
leukemia patients

The Z-residuals can be divided into k groups by cutting the linear
predictors or covariates into equally-spaced intervals. Then we can
check whether the Z-residuals of the k groups are homogeneously
distributed. A quantitative method to assess the homogeneity of such
grouped Z-residuals is to test the equality of group means or variances
of the Z-residuals. We employ the F-test in ANOVA to assess the equality
of means and Bartlett’s test to examine the equality of variances.

The scatterplots of Z-residuals against the linear predictor don’t
exhibit visible trends; their LOWESS lines are very close to the
horizontal line at 0; the boxplots of Z-residuals grouped by cutting
linear predictors into equal-spaced intervals appear to have equal means
and variance across groups. The Z-AOV and Z-BL for linear predictors
tests also gives large p-values for the wbc and lwbc models, where Z-AOV
and Z-BL are the methods of applying ANOVA and Bartlett to test the
equality of the means and variances of Z-residuals against the groups
formed with the linear predictor.

Code

``` r
for (i in 1:10) {
par(mfrow = c(2, 2), mar = c(4, 4, 1.5, 2))
plot(
    Zresid.LeukSurv.wbc,x_axis_var="index",
    main.title = "Z-residual Scatterplot of wbc model",
    irep=i
  )
plot(
    Zresid.LeukSurv.logwbc,x_axis_var="index",
    main.title = "Z-residual Scatterplot of lwbc model",
    irep=i
  )

boxplot(
    Zresid.LeukSurv.wbc,x_axis_var = "lp",
    main.title = "Z-residual Boxplot of wbc model",
    irep=i
  )
boxplot(
    Zresid.LeukSurv.logwbc,x_axis_var = "lp",
    main.title = "Z-residual Boxplot of lwbc model",
    irep=i
  )
}
```

![](zresidual_demo_files/figure-html/fig-zresid-scatter-box-.gif)

Figure 1: Figure 2: Scatter plots and box plots of Z-residuals against
LP for the wbc (left panels) and lwbc (right panels) models fitted to
the survival data of acute myeloid leukemia patients

Identifing specific model misspecifications using Z-residuals as
follows:

The above diagnostics results reveal no serious misspecification in
these two models. However, the inspection of the Z-residuals against the
covariate wbc/log(wbc) reveals that the functional form of the lwbc
model is likely misspecified. The scatterplots and comparative boxplots
of the Z-residuals against wbc/log(wbc) are shown below. The LOWESS
curve of the wbc model appears to align well with the horizontal line at
0 and the grouped Z-residuals of the wbc model appear to have equal
means and variances across groups. However, the diagnosis results for
the lwbc model are very different. It appears that there is a non-linear
trend in the LOWESS curve of the lwbc model and the grouped Z-residuals
appear to have different means across groups. The Z-AOV and Z-BL for
covariate wbc and log(wbc) also gives p-values for the wbc and lwbc
models as shown in the boxplots. The very small p-value of the Z-AOV for
covariate log(wbc) test for the lwbc models strongly suggests that the
log transformation of wbc is likely inappropriate for modelling the
survival time.

Code

``` r
for (i in 1:10) {
  par(mfrow = c(2, 2), mar = c(4, 4, 1.5, 2))

  plot(
    Zresid.LeukSurv.wbc,
    x_axis_var = "wbc",
    main.title = "Z-residual Scatterplot of wbc model",
    irep=i
  )
  plot(
    Zresid.LeukSurv.logwbc,
    x_axis_var = "logwbc",
    main.title = "Z-residual Scatterplot of lwbc model",
    irep=i
  )

  boxplot(
    Zresid.LeukSurv.wbc,
    x_axis_var = "wbc",
    main.title = "Z-residual Boxplot of wbc model",
    irep=i
  )
  boxplot(
    Zresid.LeukSurv.logwbc,
    x_axis_var = "logwbc",
    main.title = "Z-residual Boxplot of lwbc model",
    irep=i
  )
}
```

![](zresidual_demo_files/figure-html/fig-zresid-scatter-box-wbc-.gif)

Figure 2: Figure 3: Scatter plots and box plots of Z-residuals against
covariate (wbc) for the wbc (left panels) and lwbc (right panels) models
fitted to the survival data of acute myeloid leukemia patients

The boxplots of the Z-residuals against categorical covariate sex shows
the grouped Z-residuals appear to have equal means and variances across
groups. The p-values of Z-AOV and Z-BL are greater than 0.05.

### 4.5 Diagnostic Tests with Z-residuals

The Shapiro-Wilk (SW) or Shapiro-Francia (SF) normality tests applied to
Z-residuals can be used to numerically test the overall GOF of the
model. Moreover, the Shapiro-Franciard (SF) test can be employed to
assess the normality of censored data. The diagnosis of the GOF of
survival probabilities can be converted into a diagnosis of the
normality of the censored Z-residuals. Thus, by utilizing the
gofTestCensored function from the R package EnvStats, one can examine
the normality of censored Z-residuals.

Code

``` r
sw.wbc<-sw.test.zresid(Zresid.LeukSurv.wbc)
sw.lwbc<-sw.test.zresid(Zresid.LeukSurv.logwbc)
sf.wbc<-sf.test.zresid(Zresid.LeukSurv.wbc)
sf.lwbc<-sf.test.zresid(Zresid.LeukSurv.logwbc)
gof_tests<-data.frame(sw.wbc,sw.lwbc,sf.wbc,sf.lwbc)
```

The Z-residuals can be divided into k groups by cutting the covariates
or linear predictors into equally-spaced intervals. To quantitatively
evaluate the homogeneity of grouped Z-residuals, we propose testing the
equality of group means and group variances. For this purpose, we employ
the F-test in ANOVA to assess the equality of means and Bartlett’s test
to examine the equality of variances.

Code

``` r
library(Zresidual)

packageVersion("Zresidual")
```

    [1] '0.1.0'

Code

``` r
ls("package:Zresidual")[grep("bartlett", ls("package:Zresidual"))]
```

    [1] "bartlett.test.zresid"

Code

``` r
aov.wbc.lp<-aov.test.zresid(Zresid.LeukSurv.wbc,X = "lp", k.anova=10)
aov.lwbc.lp<-aov.test.zresid(Zresid.LeukSurv.logwbc,X = "lp", k.anova=10)

bl.wbc.lp<-bartlett.test.zresid(Zresid.LeukSurv.wbc,X = "lp", k.bl=10)
bl.lwbc.lp<-bartlett.test.zresid(Zresid.LeukSurv.logwbc,X = "lp", k.bl=10)

aov.wbc<-aov.test.zresid(Zresid.LeukSurv.wbc,X="wbc", k.anova=10)
aov.lwbc<-aov.test.zresid(Zresid.LeukSurv.logwbc,X="logwbc", k.anova=10)

bl.wbc<-bartlett.test.zresid(Zresid.LeukSurv.wbc,X="wbc", k.bl=10)
bl.lwbc<-bartlett.test.zresid(Zresid.LeukSurv.logwbc,X="logwbc", k.bl=10)

homogeneity_tests<-data.frame(aov.wbc.lp,aov.lwbc.lp,bl.wbc.lp,bl.lwbc.lp,
                              aov.wbc,aov.lwbc,bl.wbc,bl.lwbc)

homogeneity_tests %>%
  gt() %>%
  tab_header(
    title = "Summary of Residual Homogeneity Tests"
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 4 
  )
```

| Summary of Residual Homogeneity Tests |  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|----|
| aov.wbc.lp | aov.lwbc.lp | bl.wbc.lp | bl.lwbc.lp | aov.wbc | aov.lwbc | bl.wbc | bl.lwbc |
| 0.8701 | 0.7299 | 0.5704 | 0.0612 | 0.5868 | 0.0001 | 0.9885 | 0.2790 |
| 0.8412 | 0.9480 | 0.6401 | 0.1581 | 0.6530 | 0.0000 | 0.9172 | 0.2637 |
| 0.9724 | 0.7668 | 0.6459 | 0.7844 | 0.5199 | 0.0000 | 0.7273 | 0.5765 |
| 0.7362 | 0.8222 | 0.3540 | 0.8142 | 0.6979 | 0.0000 | 0.5464 | 0.6532 |
| 0.9054 | 0.6409 | 0.8773 | 0.2562 | 0.6199 | 0.0000 | 0.9712 | 0.1497 |
| 0.9276 | 0.7744 | 0.9451 | 0.7756 | 0.8013 | 0.0000 | 0.2140 | 0.5403 |
| 0.9730 | 0.8900 | 0.1287 | 0.1508 | 0.5196 | 0.0001 | 0.4639 | 0.9744 |
| 0.9149 | 0.8252 | 0.9055 | 0.3982 | 0.5772 | 0.0001 | 0.9324 | 0.6336 |
| 0.8917 | 0.5897 | 0.8742 | 0.1209 | 0.5520 | 0.0000 | 0.7991 | 0.9185 |
| 0.9076 | 0.6479 | 0.7593 | 0.2998 | 0.6597 | 0.0000 | 0.9469 | 0.7640 |
| 0.8491 | 0.4677 | 0.5974 | 0.7699 | 0.6240 | 0.0001 | 0.9164 | 0.8501 |
| 0.8344 | 0.6106 | 0.6324 | 0.4197 | 0.6261 | 0.0000 | 0.9602 | 0.2969 |
| 0.5545 | 0.6589 | 0.3388 | 0.5036 | 0.5046 | 0.0000 | 0.9149 | 0.1053 |
| 0.8113 | 0.6326 | 0.5928 | 0.1707 | 0.4950 | 0.0000 | 0.9720 | 0.3983 |
| 0.9349 | 0.7169 | 0.7380 | 0.4246 | 0.6726 | 0.0000 | 0.6333 | 0.2057 |
| 0.8288 | 0.4692 | 0.3657 | 0.6189 | 0.6999 | 0.0000 | 0.6881 | 0.9171 |
| 0.8953 | 0.7249 | 0.6472 | 0.7169 | 0.5148 | 0.0003 | 0.9233 | 0.4998 |
| 0.8283 | 0.6976 | 0.4643 | 0.3115 | 0.6807 | 0.0000 | 0.6673 | 0.1563 |
| 0.9857 | 0.9255 | 0.7845 | 0.7009 | 0.7497 | 0.0000 | 0.8501 | 0.2109 |
| 0.9680 | 0.7250 | 0.8388 | 0.6268 | 0.6555 | 0.0000 | 0.7037 | 0.4604 |
| 0.8550 | 0.7462 | 0.6562 | 0.2820 | 0.6430 | 0.0002 | 0.9619 | 0.7103 |
| 0.8516 | 0.4751 | 0.8086 | 0.7125 | 0.4579 | 0.0000 | 0.9993 | 0.9924 |
| 0.8516 | 0.5633 | 0.6237 | 0.9632 | 0.4680 | 0.0000 | 0.8971 | 0.4090 |
| 0.9731 | 0.7464 | 0.8468 | 0.9899 | 0.6897 | 0.0001 | 0.9002 | 0.1541 |
| 0.7465 | 0.7750 | 0.5168 | 0.2524 | 0.4849 | 0.0000 | 0.9558 | 0.6049 |
| 0.9553 | 0.8733 | 0.5593 | 0.2134 | 0.4897 | 0.0000 | 0.7720 | 0.3603 |
| 0.8583 | 0.8511 | 0.3111 | 0.9137 | 0.6133 | 0.0001 | 0.7821 | 0.0670 |
| 0.7399 | 0.9090 | 0.6836 | 0.1699 | 0.6621 | 0.0000 | 0.7433 | 0.0612 |
| 0.7602 | 0.5814 | 0.7906 | 0.1780 | 0.6619 | 0.0000 | 0.6716 | 0.3159 |
| 0.6371 | 0.8146 | 0.3231 | 0.5818 | 0.5426 | 0.0001 | 0.8762 | 0.0243 |
| 0.8832 | 0.9163 | 0.5552 | 0.1253 | 0.5224 | 0.0000 | 0.9314 | 0.9840 |
| 0.4468 | 0.4952 | 0.4761 | 0.4182 | 0.6835 | 0.0000 | 0.4865 | 0.9667 |
| 0.8672 | 0.6067 | 0.8115 | 0.8337 | 0.6370 | 0.0000 | 0.9194 | 0.3249 |
| 0.8865 | 0.9821 | 0.7892 | 0.4786 | 0.5307 | 0.0000 | 0.8697 | 0.1518 |
| 0.9592 | 0.8750 | 0.6790 | 0.5507 | 0.5010 | 0.0000 | 0.9339 | 0.4512 |
| 0.8395 | 0.7865 | 0.9027 | 0.5956 | 0.5153 | 0.0000 | 0.9601 | 0.5214 |
| 0.9496 | 0.8728 | 0.7647 | 0.8492 | 0.4693 | 0.0001 | 0.6106 | 0.5474 |
| 0.8224 | 0.8019 | 0.5463 | 0.9537 | 0.5501 | 0.0000 | 0.9517 | 0.5742 |
| 0.8000 | 0.9013 | 0.5169 | 0.6056 | 0.7378 | 0.0001 | 0.8999 | 0.9020 |
| 0.9653 | 0.8325 | 0.3091 | 0.3998 | 0.6794 | 0.0002 | 0.8375 | 0.1250 |
| 0.6236 | 0.9580 | 0.4378 | 0.4076 | 0.5674 | 0.0000 | 0.8915 | 0.1149 |
| 0.9121 | 0.7731 | 0.8035 | 0.1234 | 0.7242 | 0.0000 | 0.7681 | 0.3321 |
| 0.8218 | 0.6521 | 0.7885 | 0.1787 | 0.3675 | 0.0001 | 0.8880 | 0.9773 |
| 0.9354 | 0.8364 | 0.8702 | 0.2914 | 0.5287 | 0.0001 | 0.9823 | 0.5152 |
| 0.8635 | 0.9000 | 0.7201 | 0.4599 | 0.7043 | 0.0000 | 0.8815 | 0.4775 |
| 0.7880 | 0.8818 | 0.3463 | 0.5771 | 0.5076 | 0.0006 | 0.9578 | 0.0835 |
| 0.9727 | 0.3961 | 0.5393 | 0.4894 | 0.6594 | 0.0000 | 0.5520 | 0.1812 |
| 0.8649 | 0.8889 | 0.7808 | 0.3022 | 0.6127 | 0.0001 | 0.8477 | 0.6658 |
| 0.7780 | 0.7494 | 0.4101 | 0.3127 | 0.6202 | 0.0000 | 0.6973 | 0.9024 |
| 0.6457 | 0.7561 | 0.5731 | 0.4236 | 0.5062 | 0.0000 | 0.9956 | 0.9990 |
| 0.9193 | 0.6846 | 0.6955 | 0.9528 | 0.7055 | 0.0000 | 0.8369 | 0.3024 |
| 0.9143 | 0.8639 | 0.8141 | 0.5068 | 0.6489 | 0.0001 | 0.7016 | 0.5753 |
| 0.9102 | 0.6206 | 0.7458 | 0.5119 | 0.6184 | 0.0004 | 0.8832 | 0.0605 |
| 0.8016 | 0.8355 | 0.6097 | 0.7666 | 0.7166 | 0.0001 | 0.7388 | 0.7200 |
| 0.5824 | 0.8940 | 0.3874 | 0.4951 | 0.5847 | 0.0000 | 0.9952 | 0.5523 |
| 0.9717 | 0.7582 | 0.7118 | 0.8707 | 0.7824 | 0.0003 | 0.7240 | 0.5455 |
| 0.8691 | 0.9393 | 0.7834 | 0.4855 | 0.6889 | 0.0005 | 0.7149 | 0.2193 |
| 0.6424 | 0.6748 | 0.4293 | 0.7177 | 0.7185 | 0.0000 | 0.8373 | 0.9150 |
| 0.8288 | 0.5155 | 0.8132 | 0.2939 | 0.7340 | 0.0003 | 0.6281 | 0.4281 |
| 0.9363 | 0.5937 | 0.6353 | 0.7910 | 0.5407 | 0.0000 | 0.9798 | 0.5730 |
| 0.9697 | 0.7802 | 0.7687 | 0.1605 | 0.7100 | 0.0000 | 0.8833 | 0.6342 |
| 0.6918 | 0.8849 | 0.5202 | 0.2242 | 0.4479 | 0.0000 | 0.9144 | 0.3144 |
| 0.8829 | 0.4849 | 0.4784 | 0.7002 | 0.5697 | 0.0000 | 0.8849 | 0.7586 |
| 0.9430 | 0.7650 | 0.7832 | 0.8630 | 0.5912 | 0.0000 | 0.5266 | 0.2978 |
| 0.6914 | 0.8824 | 0.8094 | 0.1490 | 0.4066 | 0.0002 | 0.7621 | 0.0734 |
| 0.7674 | 0.8859 | 0.2502 | 0.0766 | 0.7229 | 0.0001 | 0.7342 | 0.1687 |
| 0.7079 | 0.8095 | 0.8047 | 0.4122 | 0.5776 | 0.0000 | 0.7608 | 0.7684 |
| 0.8961 | 0.5717 | 0.8127 | 0.7960 | 0.5035 | 0.0000 | 0.6655 | 0.6444 |
| 0.8625 | 0.8612 | 0.7162 | 0.6288 | 0.6473 | 0.0002 | 0.8302 | 0.6302 |
| 0.6400 | 0.5546 | 0.8139 | 0.2830 | 0.5941 | 0.0000 | 0.9286 | 0.9164 |
| 0.5144 | 0.8541 | 0.1880 | 0.4292 | 0.7357 | 0.0000 | 0.7836 | 0.8561 |
| 0.8687 | 0.7037 | 0.8541 | 0.6692 | 0.7323 | 0.0001 | 0.8757 | 0.7144 |
| 0.8853 | 0.6034 | 0.6477 | 0.2908 | 0.6283 | 0.0000 | 0.9625 | 0.8834 |
| 0.9010 | 0.8901 | 0.9007 | 0.6682 | 0.4980 | 0.0001 | 0.9584 | 0.5195 |
| 0.9155 | 0.7615 | 0.8114 | 0.5080 | 0.5753 | 0.0002 | 0.9967 | 0.8371 |
| 0.8501 | 0.9323 | 0.3364 | 0.2248 | 0.6742 | 0.0000 | 0.6223 | 0.9811 |
| 0.8731 | 0.9448 | 0.8895 | 0.1081 | 0.4279 | 0.0001 | 0.9002 | 0.1943 |
| 0.7202 | 0.3762 | 0.4979 | 0.3668 | 0.7071 | 0.0000 | 0.8072 | 0.7739 |
| 0.9043 | 0.8771 | 0.4830 | 0.2676 | 0.6775 | 0.0010 | 0.8202 | 0.1823 |
| 0.9191 | 0.8835 | 0.7672 | 0.2577 | 0.6744 | 0.0002 | 0.9815 | 0.1242 |
| 0.9635 | 0.9453 | 0.3771 | 0.1720 | 0.5728 | 0.0000 | 0.5237 | 0.7327 |
| 0.9200 | 0.9224 | 0.8046 | 0.4742 | 0.6505 | 0.0000 | 0.8645 | 0.5205 |
| 0.9362 | 0.5742 | 0.6119 | 0.6036 | 0.6099 | 0.0000 | 0.8550 | 0.5534 |
| 0.9462 | 0.7319 | 0.6994 | 0.2371 | 0.5172 | 0.0000 | 0.9875 | 0.8621 |
| 0.9354 | 0.7656 | 0.7075 | 0.3923 | 0.6707 | 0.0000 | 0.7661 | 0.0505 |
| 0.8531 | 0.7245 | 0.6293 | 0.6784 | 0.6114 | 0.0002 | 0.9512 | 0.2057 |
| 0.9302 | 0.4845 | 0.8216 | 0.3744 | 0.7238 | 0.0001 | 0.7089 | 0.5311 |
| 0.9851 | 0.8445 | 0.7617 | 0.7690 | 0.5993 | 0.0004 | 0.9675 | 0.0714 |
| 0.7232 | 0.5526 | 0.3676 | 0.2376 | 0.6238 | 0.0000 | 0.9835 | 0.3702 |
| 0.7108 | 0.7344 | 0.5674 | 0.4919 | 0.6788 | 0.0000 | 0.9066 | 0.7523 |
| 0.9528 | 0.7040 | 0.4617 | 0.6115 | 0.7855 | 0.0000 | 0.7115 | 0.8782 |
| 0.5717 | 0.7950 | 0.3309 | 0.7942 | 0.5543 | 0.0000 | 0.9678 | 0.7691 |
| 0.9263 | 0.5855 | 0.7226 | 0.9177 | 0.2596 | 0.0000 | 0.6587 | 0.0704 |
| 0.9201 | 0.6983 | 0.8365 | 0.8493 | 0.4552 | 0.0000 | 0.6933 | 0.7641 |
| 0.9529 | 0.7692 | 0.4238 | 0.3254 | 0.4073 | 0.0000 | 0.9091 | 0.8232 |
| 0.9462 | 0.6126 | 0.7488 | 0.6505 | 0.6927 | 0.0000 | 0.9403 | 0.9232 |
| 0.8367 | 0.7219 | 0.4814 | 0.6204 | 0.5948 | 0.0000 | 0.9367 | 0.7615 |
| 0.9229 | 0.5611 | 0.9221 | 0.1890 | 0.6327 | 0.0002 | 0.7034 | 0.9278 |
| 0.9017 | 0.7854 | 0.7716 | 0.6259 | 0.6064 | 0.0000 | 0.7750 | 0.4007 |
| 0.8868 | 0.7909 | 0.7685 | 0.1555 | 0.6124 | 0.0000 | 0.9301 | 0.8735 |
| 0.9654 | 0.8385 | 0.6238 | 0.8276 | 0.5886 | 0.0002 | 0.8977 | 0.2145 |
| 0.9419 | 0.8099 | 0.8102 | 0.7261 | 0.4859 | 0.0001 | 0.9519 | 0.0628 |
| 0.9161 | 0.8152 | 0.5878 | 0.4294 | 0.5187 | 0.0000 | 0.8829 | 0.2781 |
| 0.8773 | 0.6408 | 0.8637 | 0.6764 | 0.5789 | 0.0000 | 0.9104 | 0.5015 |
| 0.8803 | 0.9464 | 0.2380 | 0.5852 | 0.6107 | 0.0001 | 0.7026 | 0.1775 |
| 0.9430 | 0.8835 | 0.5206 | 0.5170 | 0.5628 | 0.0002 | 0.9844 | 0.3035 |
| 0.9161 | 0.8507 | 0.7578 | 0.1586 | 0.7019 | 0.0000 | 0.7318 | 0.7260 |
| 0.9079 | 0.5961 | 0.6708 | 0.2436 | 0.6618 | 0.0000 | 0.8592 | 0.7671 |
| 0.9023 | 0.7227 | 0.5440 | 0.6506 | 0.6621 | 0.0001 | 0.9108 | 0.3780 |
| 0.9461 | 0.7898 | 0.5727 | 0.7023 | 0.7390 | 0.0000 | 0.6131 | 0.8716 |
| 0.8862 | 0.7157 | 0.6862 | 0.6423 | 0.5666 | 0.0000 | 0.8614 | 0.0033 |
| 0.8808 | 0.8450 | 0.7311 | 0.3552 | 0.6815 | 0.0001 | 0.8789 | 0.5480 |
| 0.4912 | 0.8434 | 0.2893 | 0.2150 | 0.6199 | 0.0000 | 0.6495 | 0.5838 |
| 0.9015 | 0.7123 | 0.8693 | 0.3516 | 0.6954 | 0.0000 | 0.6867 | 0.3088 |
| 0.8036 | 0.9521 | 0.6460 | 0.2523 | 0.5756 | 0.0000 | 0.9963 | 0.9851 |
| 0.7495 | 0.9408 | 0.8052 | 0.7102 | 0.5994 | 0.0000 | 0.9375 | 0.6366 |
| 0.9313 | 0.8155 | 0.4391 | 0.2033 | 0.6525 | 0.0000 | 0.9887 | 0.1235 |
| 0.9004 | 0.6342 | 0.8044 | 0.9335 | 0.4798 | 0.0001 | 0.9349 | 0.3165 |
| 0.7133 | 0.8598 | 0.6159 | 0.7097 | 0.4751 | 0.0000 | 0.9176 | 0.8860 |
| 0.8828 | 0.7192 | 0.8307 | 0.3808 | 0.6193 | 0.0002 | 0.7077 | 0.1505 |
| 0.9369 | 0.6623 | 0.6050 | 0.4962 | 0.5357 | 0.0000 | 0.9018 | 0.2952 |
| 0.9202 | 0.7734 | 0.6529 | 0.2658 | 0.5172 | 0.0007 | 0.9983 | 0.4646 |
| 0.7536 | 0.2946 | 0.7736 | 0.5303 | 0.6686 | 0.0000 | 0.9624 | 0.8706 |
| 0.9342 | 0.9252 | 0.7708 | 0.0919 | 0.6034 | 0.0001 | 0.9724 | 0.8727 |
| 0.9156 | 0.7372 | 0.6386 | 0.8354 | 0.6257 | 0.0000 | 0.7974 | 0.9182 |
| 0.9751 | 0.7849 | 0.6965 | 0.9945 | 0.7756 | 0.0001 | 0.9139 | 0.3443 |
| 0.9218 | 0.7642 | 0.5785 | 0.3472 | 0.7482 | 0.0001 | 0.4860 | 0.0257 |
| 0.8323 | 0.7953 | 0.5383 | 0.1327 | 0.4594 | 0.0003 | 0.9852 | 0.1881 |
| 0.9315 | 0.7381 | 0.8607 | 0.1735 | 0.5863 | 0.0001 | 0.8784 | 0.1438 |
| 0.9289 | 0.5312 | 0.8465 | 0.3918 | 0.7039 | 0.0001 | 0.8044 | 0.3209 |
| 0.8882 | 0.8196 | 0.7957 | 0.6696 | 0.6809 | 0.0000 | 0.8945 | 0.3101 |
| 0.6809 | 0.8640 | 0.6429 | 0.5637 | 0.6954 | 0.0000 | 0.8381 | 0.9366 |
| 0.7274 | 0.7670 | 0.6355 | 0.4630 | 0.5059 | 0.0000 | 0.9568 | 0.7218 |
| 0.8731 | 0.6540 | 0.7860 | 0.0951 | 0.6324 | 0.0000 | 0.8646 | 0.5174 |
| 0.8643 | 0.6286 | 0.9046 | 0.4042 | 0.6303 | 0.0002 | 0.6950 | 0.8071 |
| 0.8666 | 0.6804 | 0.7965 | 0.4846 | 0.5356 | 0.0000 | 0.9749 | 0.7127 |
| 0.9149 | 0.8504 | 0.3131 | 0.7084 | 0.2409 | 0.0001 | 0.3968 | 0.5961 |
| 0.9801 | 0.8687 | 0.6828 | 0.5241 | 0.7228 | 0.0000 | 0.9022 | 0.4307 |
| 0.8181 | 0.7296 | 0.4802 | 0.6779 | 0.6759 | 0.0002 | 0.8861 | 0.5430 |
| 0.8842 | 0.6042 | 0.6835 | 0.3857 | 0.5767 | 0.0006 | 0.7946 | 0.5413 |
| 0.9641 | 0.7264 | 0.8161 | 0.2535 | 0.6369 | 0.0001 | 0.8885 | 0.2773 |
| 0.9515 | 0.7269 | 0.8751 | 0.4214 | 0.4818 | 0.0009 | 0.9384 | 0.0073 |
| 0.9043 | 0.7726 | 0.6369 | 0.7001 | 0.6339 | 0.0004 | 0.9837 | 0.1650 |
| 0.9150 | 0.4745 | 0.7323 | 0.8532 | 0.5551 | 0.0000 | 0.7430 | 0.9099 |
| 0.9501 | 0.4352 | 0.8440 | 0.2719 | 0.5910 | 0.0001 | 0.8326 | 0.7719 |
| 0.9446 | 0.6681 | 0.7762 | 0.6570 | 0.7535 | 0.0000 | 0.7517 | 0.9949 |
| 0.9944 | 0.8698 | 0.6549 | 0.1504 | 0.6931 | 0.0000 | 0.6682 | 0.2593 |
| 0.8825 | 0.6143 | 0.6966 | 0.3708 | 0.5672 | 0.0000 | 0.9483 | 0.1305 |
| 0.9099 | 0.8857 | 0.7757 | 0.2294 | 0.6806 | 0.0000 | 0.9473 | 0.2471 |
| 0.8534 | 0.7444 | 0.2961 | 0.4674 | 0.6759 | 0.0000 | 0.9573 | 0.8440 |
| 0.8998 | 0.6579 | 0.2099 | 0.8311 | 0.5673 | 0.0000 | 0.7195 | 0.0444 |
| 0.6921 | 0.9194 | 0.5631 | 0.5207 | 0.6449 | 0.0000 | 0.8634 | 0.5251 |
| 0.8838 | 0.4895 | 0.6541 | 0.4339 | 0.6561 | 0.0001 | 0.9502 | 0.9090 |
| 0.8830 | 0.8260 | 0.8584 | 0.6158 | 0.6003 | 0.0000 | 0.8653 | 0.3145 |
| 0.8283 | 0.8205 | 0.6300 | 0.1250 | 0.6863 | 0.0001 | 0.9245 | 0.0899 |
| 0.8171 | 0.8463 | 0.5318 | 0.8182 | 0.5442 | 0.0006 | 0.8813 | 0.0284 |
| 0.8882 | 0.5414 | 0.8117 | 0.3442 | 0.6324 | 0.0000 | 0.7019 | 0.7073 |
| 0.9026 | 0.9010 | 0.5243 | 0.3415 | 0.5329 | 0.0001 | 0.9795 | 0.3404 |
| 0.8401 | 0.8037 | 0.7557 | 0.7481 | 0.5445 | 0.0000 | 0.7902 | 0.1787 |
| 0.8742 | 0.9134 | 0.6697 | 0.3618 | 0.5554 | 0.0000 | 0.8454 | 0.7796 |
| 0.8667 | 0.7105 | 0.2304 | 0.6904 | 0.6343 | 0.0009 | 0.9904 | 0.5962 |
| 0.9429 | 0.8728 | 0.8608 | 0.5025 | 0.5714 | 0.0002 | 0.9570 | 0.0312 |
| 0.6857 | 0.8340 | 0.4220 | 0.2306 | 0.5868 | 0.0002 | 0.8212 | 0.5166 |
| 0.9236 | 0.7113 | 0.8111 | 0.2346 | 0.6380 | 0.0001 | 0.9368 | 0.6558 |
| 0.8798 | 0.9226 | 0.8340 | 0.6887 | 0.6533 | 0.0000 | 0.6889 | 0.8865 |
| 0.9385 | 0.7175 | 0.7743 | 0.3144 | 0.5748 | 0.0000 | 0.3263 | 0.3688 |
| 0.8000 | 0.5393 | 0.7325 | 0.6637 | 0.6495 | 0.0000 | 0.6642 | 0.1244 |
| 0.8812 | 0.8461 | 0.7697 | 0.3686 | 0.6467 | 0.0002 | 0.9787 | 0.1386 |
| 0.8475 | 0.7707 | 0.6555 | 0.9623 | 0.6397 | 0.0000 | 0.9374 | 0.6797 |
| 0.9113 | 0.8189 | 0.7719 | 0.7311 | 0.6325 | 0.0001 | 0.8463 | 0.0499 |
| 0.7451 | 0.7747 | 0.3399 | 0.5670 | 0.7529 | 0.0001 | 0.8579 | 0.3423 |
| 0.8595 | 0.9123 | 0.7030 | 0.6276 | 0.7116 | 0.0002 | 0.9071 | 0.6446 |
| 0.9459 | 0.9182 | 0.6770 | 0.5089 | 0.5638 | 0.0000 | 0.9122 | 0.9817 |
| 0.9038 | 0.8010 | 0.3597 | 0.6116 | 0.6445 | 0.0000 | 0.5676 | 0.4247 |
| 0.5625 | 0.4997 | 0.1021 | 0.7467 | 0.7289 | 0.0000 | 0.7882 | 0.8212 |
| 0.9580 | 0.6340 | 0.4334 | 0.5238 | 0.7444 | 0.0000 | 0.6738 | 0.8856 |
| 0.6581 | 0.8599 | 0.5907 | 0.9329 | 0.6010 | 0.0001 | 0.8358 | 0.3020 |
| 0.9565 | 0.8144 | 0.7213 | 0.3607 | 0.6050 | 0.0001 | 0.9041 | 0.7192 |
| 0.8729 | 0.8030 | 0.8163 | 0.5542 | 0.5753 | 0.0002 | 0.7128 | 0.8721 |
| 0.6848 | 0.7582 | 0.6090 | 0.0195 | 0.5586 | 0.0001 | 0.9888 | 0.5515 |
| 0.7789 | 0.9844 | 0.8018 | 0.2727 | 0.5624 | 0.0002 | 0.9941 | 0.0672 |
| 0.9103 | 0.9178 | 0.5983 | 0.0704 | 0.6836 | 0.0003 | 0.8165 | 0.7956 |
| 0.9248 | 0.9068 | 0.8758 | 0.4925 | 0.6159 | 0.0001 | 0.8380 | 0.6856 |
| 0.9110 | 0.3255 | 0.3939 | 0.8684 | 0.6855 | 0.0002 | 0.9445 | 0.2440 |
| 0.7863 | 0.8358 | 0.5628 | 0.7759 | 0.4953 | 0.0001 | 0.9216 | 0.1954 |
| 0.8124 | 0.8911 | 0.6759 | 0.2505 | 0.4134 | 0.0002 | 0.8979 | 0.6877 |
| 0.7577 | 0.3779 | 0.8158 | 0.9209 | 0.7227 | 0.0000 | 0.8262 | 0.5948 |
| 0.9211 | 0.8178 | 0.6389 | 0.7062 | 0.6991 | 0.0001 | 0.8049 | 0.9040 |
| 0.8965 | 0.7968 | 0.6008 | 0.0859 | 0.5575 | 0.0006 | 0.9667 | 0.2401 |
| 0.9486 | 0.6339 | 0.6275 | 0.9334 | 0.7192 | 0.0001 | 0.6539 | 0.3189 |
| 0.9523 | 0.6847 | 0.8052 | 0.4251 | 0.5055 | 0.0001 | 0.9913 | 0.6545 |
| 0.9613 | 0.6258 | 0.6635 | 0.3101 | 0.5790 | 0.0000 | 0.9714 | 0.6960 |
| 0.8100 | 0.6768 | 0.2133 | 0.6267 | 0.5617 | 0.0000 | 0.9698 | 0.8864 |
| 0.8256 | 0.6335 | 0.4930 | 0.5112 | 0.5569 | 0.0002 | 0.8284 | 0.0067 |
| 0.8772 | 0.7154 | 0.4756 | 0.2726 | 0.5463 | 0.0000 | 0.7446 | 0.4086 |
| 0.6517 | 0.7887 | 0.8434 | 0.9080 | 0.6476 | 0.0001 | 0.9651 | 0.4827 |
| 0.9311 | 0.8670 | 0.7076 | 0.1229 | 0.5221 | 0.0000 | 0.9839 | 0.7539 |
| 0.9689 | 0.7583 | 0.7651 | 0.8346 | 0.7929 | 0.0007 | 0.7232 | 0.6584 |
| 0.9298 | 0.5100 | 0.3940 | 0.2171 | 0.4088 | 0.0002 | 0.8666 | 0.3775 |
| 0.9806 | 0.6796 | 0.5506 | 0.8432 | 0.4310 | 0.0000 | 0.6689 | 0.9055 |
| 0.9472 | 0.8932 | 0.7273 | 0.1896 | 0.4683 | 0.0000 | 0.9928 | 0.5284 |
| 0.7849 | 0.9498 | 0.4170 | 0.6815 | 0.5361 | 0.0001 | 0.6167 | 0.6620 |
| 0.8879 | 0.7631 | 0.5886 | 0.1244 | 0.5902 | 0.0001 | 0.9603 | 0.4874 |
| 0.8125 | 0.8642 | 0.5903 | 0.5715 | 0.5568 | 0.0000 | 0.9502 | 0.7275 |
| 0.8789 | 0.3295 | 0.5220 | 0.5242 | 0.6247 | 0.0000 | 0.9672 | 0.9570 |
| 0.8289 | 0.8694 | 0.5390 | 0.3182 | 0.5122 | 0.0000 | 0.9135 | 0.9115 |
| 0.8593 | 0.8773 | 0.7657 | 0.2391 | 0.4918 | 0.0000 | 0.9464 | 0.9058 |
| 0.8833 | 0.8597 | 0.4011 | 0.4584 | 0.6218 | 0.0001 | 0.6196 | 0.4982 |
| 0.9275 | 0.7238 | 0.7467 | 0.6050 | 0.6159 | 0.0000 | 0.9510 | 0.4397 |
| 0.9073 | 0.9018 | 0.8518 | 0.1613 | 0.7305 | 0.0000 | 0.6728 | 0.9830 |
| 0.9708 | 0.8088 | 0.4859 | 0.7411 | 0.6704 | 0.0001 | 0.5341 | 0.7835 |
| 0.8942 | 0.7233 | 0.6265 | 0.3859 | 0.6350 | 0.0000 | 0.9191 | 0.6690 |
| 0.9413 | 0.7729 | 0.2612 | 0.3759 | 0.5993 | 0.0000 | 0.8997 | 0.0526 |
| 0.9436 | 0.7708 | 0.8283 | 0.6675 | 0.5865 | 0.0000 | 0.8934 | 0.4510 |
| 0.9581 | 0.6002 | 0.5602 | 0.4489 | 0.6090 | 0.0000 | 0.8419 | 0.2372 |
| 0.8653 | 0.7705 | 0.6170 | 0.7585 | 0.5723 | 0.0000 | 0.7656 | 0.7548 |
| 0.9125 | 0.9425 | 0.3969 | 0.2070 | 0.5473 | 0.0007 | 0.9012 | 0.2431 |
| 0.9314 | 0.7896 | 0.6960 | 0.7447 | 0.5942 | 0.0000 | 0.9920 | 0.3636 |
| 0.7595 | 0.7162 | 0.6525 | 0.1890 | 0.7009 | 0.0004 | 0.9777 | 0.7633 |
| 0.7456 | 0.8795 | 0.5581 | 0.3838 | 0.6508 | 0.0002 | 0.9879 | 0.7061 |
| 0.9358 | 0.6819 | 0.5025 | 0.1967 | 0.8624 | 0.0000 | 0.2987 | 0.9316 |
| 0.9101 | 0.9087 | 0.6470 | 0.5786 | 0.4470 | 0.0000 | 0.9295 | 0.3880 |
| 0.9260 | 0.8088 | 0.5565 | 0.2299 | 0.6160 | 0.0000 | 0.9699 | 0.9897 |
| 0.9012 | 0.6734 | 0.7963 | 0.6184 | 0.5469 | 0.0000 | 0.9028 | 0.7043 |
| 0.9371 | 0.8306 | 0.6775 | 0.4038 | 0.5496 | 0.0001 | 0.9595 | 0.2136 |
| 0.7643 | 0.6613 | 0.7777 | 0.4683 | 0.6046 | 0.0001 | 0.9879 | 0.6827 |
| 0.7675 | 0.8722 | 0.2568 | 0.8803 | 0.5319 | 0.0005 | 0.9547 | 0.8407 |
| 0.8103 | 0.9432 | 0.8714 | 0.6170 | 0.6782 | 0.0000 | 0.8483 | 0.9710 |
| 0.8215 | 0.6673 | 0.6963 | 0.2765 | 0.4430 | 0.0000 | 0.9842 | 0.7671 |
| 0.9153 | 0.7154 | 0.8246 | 0.7648 | 0.5074 | 0.0000 | 0.9599 | 0.8751 |
| 0.8295 | 0.5112 | 0.5540 | 0.8530 | 0.5704 | 0.0001 | 0.7917 | 0.7578 |
| 0.8231 | 0.5870 | 0.6122 | 0.3321 | 0.6171 | 0.0000 | 0.8439 | 0.8255 |
| 0.6274 | 0.5391 | 0.7368 | 0.2997 | 0.6801 | 0.0001 | 0.8431 | 0.7111 |
| 0.9493 | 0.7696 | 0.8092 | 0.2074 | 0.6970 | 0.0002 | 0.8278 | 0.3931 |
| 0.9740 | 0.6519 | 0.8609 | 0.1988 | 0.5653 | 0.0000 | 0.8454 | 0.3377 |
| 0.9220 | 0.8597 | 0.9182 | 0.0688 | 0.6864 | 0.0000 | 0.8720 | 0.8943 |
| 0.7841 | 0.8546 | 0.7495 | 0.6162 | 0.6074 | 0.0003 | 0.9675 | 0.0421 |
| 0.9634 | 0.7921 | 0.7928 | 0.6351 | 0.7186 | 0.0000 | 0.8578 | 0.3081 |
| 0.8583 | 0.8373 | 0.6976 | 0.2607 | 0.6148 | 0.0000 | 0.9874 | 0.2169 |
| 0.9137 | 0.8944 | 0.7071 | 0.3203 | 0.5878 | 0.0001 | 0.9841 | 0.1091 |
| 0.8791 | 0.8820 | 0.8135 | 0.6866 | 0.4265 | 0.0006 | 0.9145 | 0.0431 |
| 0.5147 | 0.6758 | 0.6339 | 0.6217 | 0.6616 | 0.0001 | 0.7557 | 0.1565 |
| 0.8072 | 0.6860 | 0.4845 | 0.8313 | 0.6570 | 0.0005 | 0.7492 | 0.0555 |
| 0.9659 | 0.7488 | 0.8411 | 0.8635 | 0.5570 | 0.0000 | 0.8476 | 0.5285 |
| 0.9215 | 0.8866 | 0.7053 | 0.7737 | 0.5475 | 0.0006 | 0.7430 | 0.3853 |
| 0.7939 | 0.7304 | 0.6746 | 0.4448 | 0.5813 | 0.0001 | 0.9895 | 0.4842 |
| 0.8938 | 0.6707 | 0.9322 | 0.4164 | 0.7630 | 0.0000 | 0.7869 | 0.4957 |
| 0.7315 | 0.9441 | 0.4918 | 0.0277 | 0.6599 | 0.0001 | 0.7767 | 0.3299 |
| 0.8681 | 0.7465 | 0.7702 | 0.7783 | 0.5874 | 0.0000 | 0.9690 | 0.5223 |
| 0.8622 | 0.8923 | 0.7153 | 0.7226 | 0.4945 | 0.0001 | 0.9224 | 0.8593 |
| 0.9432 | 0.8746 | 0.9300 | 0.0481 | 0.5963 | 0.0000 | 0.9580 | 0.7343 |
| 0.9723 | 0.7091 | 0.7826 | 0.8803 | 0.5189 | 0.0000 | 0.8888 | 0.3118 |
| 0.8483 | 0.9611 | 0.4760 | 0.2590 | 0.6007 | 0.0000 | 0.9099 | 0.7493 |
| 0.9163 | 0.7216 | 0.5295 | 0.8910 | 0.6304 | 0.0000 | 0.8807 | 0.2571 |
| 0.9304 | 0.8137 | 0.7141 | 0.8318 | 0.6486 | 0.0000 | 0.8323 | 0.9360 |
| 0.9399 | 0.4159 | 0.8452 | 0.9931 | 0.6345 | 0.0002 | 0.5619 | 0.6504 |
| 0.9803 | 0.5456 | 0.4625 | 0.7100 | 0.5493 | 0.0000 | 0.7849 | 0.9964 |
| 0.9364 | 0.5960 | 0.6422 | 0.7922 | 0.6431 | 0.0001 | 0.9634 | 0.5595 |
| 0.7965 | 0.3271 | 0.5517 | 0.4794 | 0.4457 | 0.0000 | 0.7606 | 0.0689 |
| 0.7376 | 0.6344 | 0.7899 | 0.0796 | 0.6359 | 0.0005 | 0.7428 | 0.0382 |
| 0.8192 | 0.6834 | 0.5971 | 0.3250 | 0.6607 | 0.0001 | 0.9798 | 0.6423 |
| 0.9269 | 0.9664 | 0.6242 | 0.0930 | 0.4751 | 0.0000 | 0.9540 | 0.6317 |
| 0.9086 | 0.5617 | 0.7468 | 0.5800 | 0.6144 | 0.0001 | 0.8791 | 0.8141 |
| 0.9708 | 0.7108 | 0.6028 | 0.5307 | 0.6039 | 0.0000 | 0.6501 | 0.3041 |
| 0.9731 | 0.5594 | 0.5436 | 0.6830 | 0.5557 | 0.0000 | 0.9803 | 0.7241 |
| 0.8373 | 0.9699 | 0.8957 | 0.4043 | 0.6301 | 0.0000 | 0.9597 | 0.8149 |
| 0.8103 | 0.6682 | 0.5524 | 0.7599 | 0.5975 | 0.0000 | 0.9730 | 0.7223 |
| 0.7979 | 0.8670 | 0.7114 | 0.7022 | 0.6424 | 0.0000 | 0.7784 | 0.7902 |
| 0.6852 | 0.8794 | 0.1587 | 0.4184 | 0.6015 | 0.0000 | 0.9889 | 0.0848 |
| 0.8633 | 0.6229 | 0.4861 | 0.3865 | 0.5967 | 0.0001 | 0.7649 | 0.6828 |
| 0.8407 | 0.6509 | 0.8431 | 0.4930 | 0.6364 | 0.0004 | 0.9111 | 0.0413 |
| 0.8793 | 0.4549 | 0.8746 | 0.2481 | 0.4370 | 0.0007 | 0.8128 | 0.7444 |
| 0.7539 | 0.9049 | 0.5873 | 0.1299 | 0.5483 | 0.0000 | 0.9684 | 0.9744 |
| 0.9180 | 0.9242 | 0.7834 | 0.6295 | 0.5998 | 0.0000 | 0.9587 | 0.0899 |
| 0.8231 | 0.5508 | 0.6504 | 0.9121 | 0.6228 | 0.0000 | 0.7885 | 0.6667 |
| 0.9213 | 0.7826 | 0.7342 | 0.2980 | 0.6247 | 0.0003 | 0.9905 | 0.2201 |
| 0.8697 | 0.7382 | 0.6934 | 0.5891 | 0.6006 | 0.0000 | 0.9179 | 0.2501 |
| 0.8909 | 0.6764 | 0.8517 | 0.4649 | 0.5287 | 0.0003 | 0.8038 | 0.6175 |
| 0.9146 | 0.3828 | 0.4871 | 0.5870 | 0.7679 | 0.0004 | 0.8212 | 0.6501 |
| 0.8700 | 0.7510 | 0.6779 | 0.8752 | 0.6243 | 0.0000 | 0.8365 | 0.4073 |
| 0.8757 | 0.9112 | 0.6511 | 0.4606 | 0.6905 | 0.0005 | 0.9006 | 0.8247 |
| 0.7372 | 0.7837 | 0.7786 | 0.4118 | 0.5464 | 0.0000 | 0.7499 | 0.3475 |
| 0.9296 | 0.9106 | 0.8547 | 0.3956 | 0.5673 | 0.0000 | 0.6874 | 0.3662 |
| 0.9254 | 0.9502 | 0.7810 | 0.2025 | 0.4876 | 0.0000 | 0.9386 | 0.4145 |
| 0.9240 | 0.9537 | 0.7320 | 0.5617 | 0.5360 | 0.0000 | 0.8525 | 0.4808 |
| 0.7375 | 0.9783 | 0.2743 | 0.0636 | 0.5898 | 0.0001 | 0.7961 | 0.4227 |
| 0.7492 | 0.9378 | 0.5947 | 0.3902 | 0.5746 | 0.0001 | 0.8981 | 0.8107 |
| 0.8811 | 0.8261 | 0.6269 | 0.4091 | 0.7870 | 0.0001 | 0.6942 | 0.3102 |
| 0.8490 | 0.6080 | 0.8142 | 0.7305 | 0.4607 | 0.0000 | 0.9852 | 0.2833 |
| 0.8205 | 0.4343 | 0.6039 | 0.3440 | 0.5923 | 0.0000 | 0.6532 | 0.8377 |
| 0.8835 | 0.9302 | 0.4966 | 0.4355 | 0.5812 | 0.0001 | 0.8666 | 0.7333 |
| 0.8719 | 0.7106 | 0.7183 | 0.2608 | 0.5779 | 0.0004 | 0.9557 | 0.0777 |
| 0.7507 | 0.4064 | 0.4845 | 0.5840 | 0.2680 | 0.0000 | 0.5562 | 0.7023 |
| 0.8942 | 0.8506 | 0.3709 | 0.4622 | 0.6747 | 0.0000 | 0.5979 | 0.9751 |
| 0.7419 | 0.5279 | 0.6910 | 0.5831 | 0.6903 | 0.0001 | 0.9164 | 0.0201 |
| 0.9304 | 0.8352 | 0.8353 | 0.5314 | 0.4957 | 0.0000 | 0.9298 | 0.3832 |
| 0.6932 | 0.8246 | 0.5428 | 0.3375 | 0.5299 | 0.0000 | 0.9895 | 0.9032 |
| 0.9439 | 0.8557 | 0.8541 | 0.6513 | 0.6266 | 0.0000 | 0.9701 | 0.9785 |
| 0.8828 | 0.7429 | 0.1839 | 0.2553 | 0.5694 | 0.0000 | 0.8906 | 0.8274 |
| 0.8699 | 0.5205 | 0.5783 | 0.1087 | 0.5166 | 0.0015 | 0.8949 | 0.2753 |
| 0.8763 | 0.6443 | 0.6531 | 0.9700 | 0.6526 | 0.0002 | 0.8202 | 0.3697 |
| 0.9353 | 0.7903 | 0.4790 | 0.4927 | 0.8024 | 0.0000 | 0.6491 | 0.4528 |
| 0.8362 | 0.5675 | 0.6435 | 0.8533 | 0.6762 | 0.0002 | 0.8188 | 0.7459 |
| 0.9643 | 0.9020 | 0.8369 | 0.5087 | 0.6084 | 0.0000 | 0.8932 | 0.5805 |
| 0.8048 | 0.8634 | 0.6600 | 0.5090 | 0.4625 | 0.0004 | 0.8782 | 0.5139 |
| 0.8137 | 0.9244 | 0.8271 | 0.4635 | 0.4589 | 0.0000 | 0.9308 | 0.7368 |
| 0.9603 | 0.8758 | 0.8351 | 0.4677 | 0.5620 | 0.0000 | 0.9080 | 0.1402 |
| 0.7669 | 0.7263 | 0.3484 | 0.4177 | 0.5840 | 0.0000 | 0.5794 | 0.6134 |
| 0.8714 | 0.8100 | 0.6152 | 0.6142 | 0.4682 | 0.0000 | 0.9643 | 0.6750 |
| 0.5241 | 0.8877 | 0.2478 | 0.0371 | 0.5573 | 0.0000 | 0.9993 | 0.9248 |
| 0.9476 | 0.8337 | 0.6604 | 0.4807 | 0.5830 | 0.0000 | 0.9857 | 0.7508 |
| 0.9148 | 0.7898 | 0.4758 | 0.0583 | 0.6711 | 0.0000 | 0.8485 | 0.4461 |
| 0.8290 | 0.8570 | 0.3025 | 0.2730 | 0.7578 | 0.0000 | 0.9017 | 0.9756 |
| 0.9572 | 0.6467 | 0.8300 | 0.5341 | 0.5461 | 0.0000 | 0.7728 | 0.9704 |
| 0.9706 | 0.8257 | 0.8620 | 0.2693 | 0.6231 | 0.0001 | 0.9582 | 0.4590 |
| 0.8236 | 0.8397 | 0.5811 | 0.1485 | 0.8140 | 0.0000 | 0.6995 | 0.3625 |
| 0.9537 | 0.5841 | 0.6468 | 0.1547 | 0.5666 | 0.0000 | 0.8242 | 0.9438 |
| 0.6153 | 0.8838 | 0.4199 | 0.2559 | 0.7300 | 0.0002 | 0.4152 | 0.1345 |
| 0.9542 | 0.8327 | 0.6365 | 0.5464 | 0.4925 | 0.0000 | 0.9664 | 0.9613 |
| 0.9307 | 0.8402 | 0.7256 | 0.5907 | 0.4999 | 0.0000 | 0.9808 | 0.6959 |
| 0.8731 | 0.9197 | 0.5081 | 0.5413 | 0.6191 | 0.0001 | 0.8214 | 0.9113 |
| 0.7856 | 0.8133 | 0.8235 | 0.9204 | 0.5409 | 0.0001 | 0.9497 | 0.2612 |
| 0.8600 | 0.9358 | 0.8019 | 0.2860 | 0.6317 | 0.0000 | 0.7694 | 0.9649 |
| 0.9322 | 0.6063 | 0.6674 | 0.2325 | 0.5950 | 0.0000 | 0.9598 | 0.8862 |
| 0.7651 | 0.4882 | 0.7903 | 0.7405 | 0.6651 | 0.0000 | 0.9528 | 0.8368 |
| 0.8908 | 0.9185 | 0.2744 | 0.9664 | 0.7013 | 0.0001 | 0.7961 | 0.3493 |
| 0.8128 | 0.8162 | 0.5070 | 0.8661 | 0.6195 | 0.0000 | 0.7839 | 0.5938 |
| 0.9098 | 0.8911 | 0.4923 | 0.4336 | 0.6213 | 0.0000 | 0.8579 | 0.9418 |
| 0.7689 | 0.6827 | 0.5457 | 0.3969 | 0.4968 | 0.0005 | 0.9934 | 0.5473 |
| 0.9281 | 0.7312 | 0.6131 | 0.7121 | 0.6330 | 0.0000 | 0.7788 | 0.1585 |
| 0.8532 | 0.6270 | 0.3869 | 0.3891 | 0.5714 | 0.0001 | 0.9758 | 0.8849 |
| 0.8356 | 0.7212 | 0.8260 | 0.7469 | 0.6461 | 0.0000 | 0.9520 | 0.9898 |
| 0.8394 | 0.6762 | 0.8493 | 0.7912 | 0.4633 | 0.0000 | 0.9600 | 0.5684 |
| 0.9179 | 0.8650 | 0.8957 | 0.4241 | 0.6969 | 0.0001 | 0.5973 | 0.4016 |
| 0.9518 | 0.5190 | 0.5129 | 0.4252 | 0.6654 | 0.0001 | 0.8227 | 0.2677 |
| 0.8293 | 0.5484 | 0.6906 | 0.7383 | 0.6317 | 0.0000 | 0.9432 | 0.6760 |
| 0.9335 | 0.8278 | 0.9394 | 0.6501 | 0.6332 | 0.0000 | 0.9751 | 0.6770 |
| 0.8369 | 0.8437 | 0.3977 | 0.1710 | 0.5910 | 0.0001 | 0.8220 | 0.9584 |
| 0.9150 | 0.5014 | 0.5833 | 0.0969 | 0.6633 | 0.0000 | 0.7156 | 0.4367 |
| 0.9283 | 0.9019 | 0.8199 | 0.4755 | 0.5956 | 0.0000 | 0.9921 | 0.9429 |
| 0.9489 | 0.7317 | 0.5852 | 0.7153 | 0.6154 | 0.0000 | 0.9259 | 0.6356 |
| 0.9652 | 0.7351 | 0.7949 | 0.0140 | 0.6750 | 0.0001 | 0.9062 | 0.1102 |
| 0.9663 | 0.8010 | 0.8847 | 0.8310 | 0.6505 | 0.0001 | 0.9104 | 0.7233 |
| 0.8445 | 0.8222 | 0.3338 | 0.5614 | 0.6645 | 0.0000 | 0.6181 | 0.7063 |
| 0.9172 | 0.7556 | 0.9102 | 0.5887 | 0.7006 | 0.0000 | 0.8048 | 0.3131 |
| 0.9555 | 0.8641 | 0.5743 | 0.6392 | 0.4677 | 0.0002 | 0.9153 | 0.6381 |
| 0.7451 | 0.5391 | 0.3915 | 0.2416 | 0.5747 | 0.0000 | 0.9307 | 0.9928 |
| 0.9034 | 0.7197 | 0.6381 | 0.2922 | 0.7973 | 0.0000 | 0.5019 | 0.4566 |
| 0.9146 | 0.6395 | 0.6967 | 0.2974 | 0.4939 | 0.0000 | 0.9021 | 0.9310 |
| 0.6937 | 0.8441 | 0.6476 | 0.1783 | 0.3199 | 0.0000 | 0.8682 | 0.9070 |
| 0.9346 | 0.9203 | 0.9060 | 0.6483 | 0.6260 | 0.0000 | 0.7679 | 0.2922 |
| 0.7743 | 0.8310 | 0.6924 | 0.2420 | 0.5287 | 0.0001 | 0.9740 | 0.9066 |
| 0.8107 | 0.9810 | 0.5164 | 0.1405 | 0.5503 | 0.0000 | 0.8180 | 0.1278 |
| 0.9131 | 0.5947 | 0.7802 | 0.4586 | 0.6629 | 0.0000 | 0.9826 | 0.7982 |
| 0.9076 | 0.9005 | 0.6113 | 0.6617 | 0.5584 | 0.0001 | 0.7999 | 0.9140 |
| 0.9365 | 0.7377 | 0.8037 | 0.1738 | 0.6739 | 0.0000 | 0.6223 | 0.2700 |
| 0.9267 | 0.7960 | 0.4257 | 0.2492 | 0.5965 | 0.0001 | 0.9340 | 0.1202 |
| 0.9223 | 0.7872 | 0.6795 | 0.2382 | 0.5036 | 0.0000 | 0.8917 | 0.0604 |
| 0.8322 | 0.7430 | 0.6787 | 0.8158 | 0.6404 | 0.0000 | 0.9394 | 0.9703 |
| 0.8086 | 0.7227 | 0.4943 | 0.9133 | 0.6488 | 0.0000 | 0.9473 | 0.6857 |
| 0.9444 | 0.5771 | 0.6658 | 0.7070 | 0.6333 | 0.0010 | 0.8200 | 0.4868 |
| 0.8069 | 0.5507 | 0.7169 | 0.2355 | 0.5112 | 0.0000 | 0.8825 | 0.9597 |
| 0.9251 | 0.5585 | 0.2523 | 0.9521 | 0.7604 | 0.0002 | 0.7157 | 0.7361 |
| 0.6701 | 0.7798 | 0.4387 | 0.7580 | 0.5937 | 0.0002 | 0.9193 | 0.6050 |
| 0.9599 | 0.7206 | 0.8933 | 0.7008 | 0.6146 | 0.0001 | 0.9696 | 0.5516 |
| 0.7703 | 0.6804 | 0.7185 | 0.6817 | 0.5545 | 0.0000 | 0.8845 | 0.9892 |
| 0.8131 | 0.5751 | 0.4850 | 0.3701 | 0.7500 | 0.0000 | 0.9197 | 0.9042 |
| 0.9757 | 0.8044 | 0.3968 | 0.7321 | 0.7031 | 0.0001 | 0.8887 | 0.5384 |
| 0.9038 | 0.8311 | 0.7464 | 0.1991 | 0.5800 | 0.0000 | 0.6683 | 0.8960 |
| 0.9589 | 0.5343 | 0.4955 | 0.4655 | 0.6343 | 0.0004 | 0.6225 | 0.3800 |
| 0.9166 | 0.8009 | 0.6044 | 0.2482 | 0.6209 | 0.0000 | 0.8418 | 0.4548 |
| 0.7658 | 0.9165 | 0.8715 | 0.2764 | 0.4608 | 0.0000 | 0.8783 | 0.4864 |
| 0.9663 | 0.9410 | 0.5124 | 0.2961 | 0.6291 | 0.0001 | 0.9223 | 0.8604 |
| 0.8685 | 0.7400 | 0.8474 | 0.6419 | 0.6904 | 0.0001 | 0.7086 | 0.8218 |
| 0.8343 | 0.5333 | 0.6383 | 0.9218 | 0.5061 | 0.0000 | 0.9579 | 0.0709 |
| 0.8857 | 0.5319 | 0.8502 | 0.8053 | 0.6186 | 0.0001 | 0.9776 | 0.4431 |
| 0.8972 | 0.5405 | 0.8731 | 0.5246 | 0.7701 | 0.0000 | 0.6923 | 0.5569 |
| 0.9555 | 0.2773 | 0.7271 | 0.0749 | 0.7301 | 0.0000 | 0.5210 | 0.9971 |
| 0.6567 | 0.8835 | 0.2890 | 0.7680 | 0.4802 | 0.0006 | 0.6826 | 0.1814 |
| 0.8901 | 0.6907 | 0.3425 | 0.6925 | 0.4683 | 0.0002 | 0.8782 | 0.7694 |
| 0.9283 | 0.7322 | 0.7454 | 0.5730 | 0.6657 | 0.0000 | 0.8592 | 0.8759 |
| 0.9766 | 0.7418 | 0.5873 | 0.8325 | 0.5685 | 0.0000 | 0.7128 | 0.9427 |
| 0.8986 | 0.7663 | 0.7163 | 0.2578 | 0.4767 | 0.0000 | 0.9810 | 0.9700 |
| 0.7517 | 0.6152 | 0.3967 | 0.6975 | 0.7202 | 0.0000 | 0.7465 | 0.6268 |
| 0.9117 | 0.5642 | 0.7303 | 0.9089 | 0.7203 | 0.0000 | 0.9365 | 0.0964 |
| 0.8621 | 0.8614 | 0.8026 | 0.5980 | 0.4731 | 0.0000 | 0.8847 | 0.4572 |
| 0.8689 | 0.7280 | 0.7194 | 0.2253 | 0.7095 | 0.0000 | 0.8385 | 0.8544 |
| 0.9131 | 0.6535 | 0.8613 | 0.9468 | 0.7225 | 0.0000 | 0.6277 | 0.7070 |
| 0.8984 | 0.6058 | 0.6572 | 0.7029 | 0.7102 | 0.0001 | 0.8761 | 0.6640 |
| 0.9114 | 0.8898 | 0.8131 | 0.6341 | 0.4938 | 0.0000 | 0.8923 | 0.3473 |
| 0.8933 | 0.8450 | 0.7338 | 0.3863 | 0.6116 | 0.0002 | 0.7955 | 0.6567 |
| 0.8514 | 0.5505 | 0.5918 | 0.6435 | 0.5530 | 0.0005 | 0.9923 | 0.7277 |
| 0.7736 | 0.7291 | 0.3172 | 0.5198 | 0.5360 | 0.0000 | 0.8669 | 0.8005 |
| 0.9611 | 0.6010 | 0.5456 | 0.4314 | 0.7080 | 0.0000 | 0.7836 | 0.4095 |
| 0.9528 | 0.8350 | 0.6698 | 0.4854 | 0.7269 | 0.0000 | 0.7600 | 0.7649 |
| 0.8494 | 0.6896 | 0.5863 | 0.3307 | 0.6351 | 0.0000 | 0.5852 | 0.7690 |
| 0.7332 | 0.7717 | 0.7705 | 0.6613 | 0.7178 | 0.0000 | 0.6132 | 0.8799 |
| 0.9057 | 0.9142 | 0.8071 | 0.6745 | 0.5521 | 0.0000 | 0.9849 | 0.9530 |
| 0.9854 | 0.4210 | 0.8079 | 0.9356 | 0.6697 | 0.0000 | 0.9341 | 0.1544 |
| 0.7491 | 0.9067 | 0.5924 | 0.2819 | 0.5870 | 0.0000 | 0.9844 | 0.5958 |
| 0.8720 | 0.7158 | 0.6887 | 0.9517 | 0.5511 | 0.0001 | 0.9772 | 0.6516 |
| 0.8470 | 0.8810 | 0.7879 | 0.2359 | 0.5467 | 0.0001 | 0.9830 | 0.0823 |
| 0.8860 | 0.8187 | 0.2797 | 0.7011 | 0.7368 | 0.0000 | 0.8427 | 0.6239 |
| 0.8998 | 0.8016 | 0.6498 | 0.6916 | 0.7366 | 0.0002 | 0.4377 | 0.3185 |
| 0.9419 | 0.7657 | 0.1582 | 0.6137 | 0.6270 | 0.0000 | 0.8407 | 0.0139 |
| 0.8724 | 0.9186 | 0.8332 | 0.1833 | 0.5097 | 0.0001 | 0.7130 | 0.5626 |
| 0.8478 | 0.8758 | 0.4852 | 0.0654 | 0.6089 | 0.0001 | 0.9693 | 0.4086 |
| 0.7080 | 0.7820 | 0.4436 | 0.5357 | 0.6409 | 0.0001 | 0.8744 | 0.2058 |
| 0.9867 | 0.8059 | 0.7306 | 0.7147 | 0.4448 | 0.0000 | 0.8341 | 0.3133 |
| 0.9467 | 0.7651 | 0.7774 | 0.6312 | 0.4921 | 0.0000 | 0.8759 | 0.1598 |
| 0.9361 | 0.9047 | 0.8989 | 0.3658 | 0.6423 | 0.0001 | 0.8843 | 0.2215 |
| 0.7301 | 0.8141 | 0.2725 | 0.1444 | 0.5518 | 0.0000 | 0.9208 | 0.3310 |
| 0.6974 | 0.7639 | 0.2123 | 0.3882 | 0.5960 | 0.0001 | 0.8409 | 0.1269 |
| 0.8005 | 0.6671 | 0.8450 | 0.7095 | 0.5582 | 0.0000 | 0.7517 | 0.6270 |
| 0.8077 | 0.9818 | 0.7329 | 0.0869 | 0.5391 | 0.0000 | 0.9069 | 0.7219 |
| 0.8934 | 0.6189 | 0.7306 | 0.8650 | 0.5136 | 0.0002 | 0.8597 | 0.3708 |
| 0.9341 | 0.9505 | 0.1969 | 0.5351 | 0.6388 | 0.0000 | 0.9433 | 0.8812 |
| 0.7371 | 0.8282 | 0.3747 | 0.7083 | 0.6308 | 0.0003 | 0.9601 | 0.5717 |
| 0.9367 | 0.9156 | 0.7445 | 0.6384 | 0.7679 | 0.0001 | 0.8932 | 0.1639 |
| 0.7645 | 0.7840 | 0.4688 | 0.3377 | 0.6313 | 0.0000 | 0.8967 | 0.6225 |
| 0.8336 | 0.7155 | 0.8322 | 0.7385 | 0.6695 | 0.0008 | 0.7143 | 0.1749 |
| 0.9362 | 0.5231 | 0.2681 | 0.8670 | 0.7009 | 0.0000 | 0.9143 | 0.4680 |
| 0.9596 | 0.7979 | 0.2900 | 0.7807 | 0.6779 | 0.0002 | 0.9410 | 0.5911 |
| 0.4650 | 0.9721 | 0.2676 | 0.1380 | 0.4790 | 0.0002 | 0.9160 | 0.1441 |
| 0.7995 | 0.6732 | 0.7708 | 0.3974 | 0.5061 | 0.0000 | 0.9853 | 0.2658 |
| 0.8979 | 0.9366 | 0.8405 | 0.5270 | 0.7107 | 0.0000 | 0.7810 | 0.8953 |
| 0.8234 | 0.8086 | 0.7809 | 0.4469 | 0.5986 | 0.0000 | 0.6576 | 0.8445 |
| 0.8186 | 0.3485 | 0.6928 | 0.4575 | 0.5999 | 0.0000 | 0.9175 | 0.4446 |
| 0.8931 | 0.7387 | 0.4192 | 0.2461 | 0.7580 | 0.0002 | 0.5812 | 0.5223 |
| 0.7664 | 0.8381 | 0.5932 | 0.3352 | 0.5119 | 0.0002 | 0.9143 | 0.3675 |
| 0.9744 | 0.5968 | 0.5963 | 0.5072 | 0.6569 | 0.0000 | 0.9592 | 0.9862 |
| 0.8075 | 0.6753 | 0.5513 | 0.3526 | 0.6118 | 0.0000 | 0.9782 | 0.9940 |
| 0.8890 | 0.4958 | 0.6392 | 0.2040 | 0.6050 | 0.0000 | 0.8015 | 0.5217 |
| 0.9647 | 0.7334 | 0.8989 | 0.8575 | 0.5948 | 0.0000 | 0.9945 | 0.7807 |
| 0.9610 | 0.6686 | 0.6877 | 0.6295 | 0.5533 | 0.0000 | 0.9719 | 0.8377 |
| 0.9963 | 0.6395 | 0.1988 | 0.3548 | 0.8484 | 0.0000 | 0.4398 | 0.2210 |
| 0.6546 | 0.5941 | 0.2403 | 0.6724 | 0.5737 | 0.0003 | 0.6997 | 0.7188 |
| 0.6804 | 0.7302 | 0.5777 | 0.5372 | 0.6621 | 0.0000 | 0.9786 | 0.9628 |
| 0.5810 | 0.8672 | 0.5007 | 0.6496 | 0.5461 | 0.0000 | 0.9998 | 0.1646 |
| 0.6859 | 0.8627 | 0.5777 | 0.1219 | 0.6109 | 0.0004 | 0.6926 | 0.3552 |
| 0.9357 | 0.8791 | 0.7195 | 0.3122 | 0.5628 | 0.0000 | 0.5021 | 0.6535 |
| 0.9313 | 0.9242 | 0.8768 | 0.1025 | 0.5498 | 0.0001 | 0.9888 | 0.2435 |
| 0.7440 | 0.4316 | 0.0732 | 0.8845 | 0.6788 | 0.0000 | 0.9544 | 0.4471 |
| 0.8831 | 0.8491 | 0.7378 | 0.0532 | 0.6805 | 0.0000 | 0.6641 | 0.9722 |
| 0.8670 | 0.4159 | 0.5709 | 0.6000 | 0.5352 | 0.0000 | 0.9672 | 0.1065 |
| 0.9363 | 0.9669 | 0.6454 | 0.6925 | 0.6773 | 0.0011 | 0.8665 | 0.0319 |
| 0.6086 | 0.7410 | 0.3162 | 0.2283 | 0.6651 | 0.0000 | 0.6521 | 0.8465 |
| 0.9638 | 0.5175 | 0.7059 | 0.8022 | 0.5075 | 0.0000 | 0.9002 | 0.7922 |
| 0.9020 | 0.5361 | 0.7367 | 0.7113 | 0.6925 | 0.0000 | 0.6980 | 0.9594 |
| 0.8935 | 0.7258 | 0.7317 | 0.6416 | 0.6268 | 0.0000 | 0.8857 | 0.5301 |
| 0.8759 | 0.7336 | 0.6874 | 0.8675 | 0.5439 | 0.0000 | 0.8391 | 0.8470 |
| 0.9269 | 0.9143 | 0.8604 | 0.2072 | 0.6661 | 0.0002 | 0.9501 | 0.1241 |
| 0.8552 | 0.6410 | 0.7311 | 0.6016 | 0.6788 | 0.0000 | 0.8245 | 0.9949 |
| 0.7195 | 0.7966 | 0.5636 | 0.6782 | 0.7707 | 0.0001 | 0.8582 | 0.3495 |
| 0.9707 | 0.7509 | 0.5824 | 0.8154 | 0.6702 | 0.0004 | 0.7842 | 0.5497 |
| 0.9587 | 0.9044 | 0.8551 | 0.4250 | 0.6203 | 0.0000 | 0.9763 | 0.1299 |
| 0.8980 | 0.8162 | 0.5690 | 0.1822 | 0.6794 | 0.0003 | 0.8610 | 0.3253 |
| 0.8256 | 0.7287 | 0.4702 | 0.6266 | 0.5561 | 0.0000 | 0.9122 | 0.6536 |
| 0.9477 | 0.9329 | 0.5315 | 0.8308 | 0.6198 | 0.0000 | 0.9193 | 0.8248 |
| 0.9518 | 0.8630 | 0.7765 | 0.7754 | 0.6679 | 0.0001 | 0.8874 | 0.0857 |
| 0.9777 | 0.8165 | 0.7144 | 0.4514 | 0.6662 | 0.0003 | 0.9740 | 0.0060 |
| 0.9130 | 0.6870 | 0.7955 | 0.2635 | 0.6405 | 0.0000 | 0.9214 | 0.3705 |
| 0.8122 | 0.5871 | 0.8365 | 0.5075 | 0.6265 | 0.0000 | 0.9974 | 0.8852 |
| 0.9780 | 0.5897 | 0.8142 | 0.7357 | 0.6239 | 0.0000 | 0.9870 | 0.4811 |
| 0.9330 | 0.6654 | 0.2579 | 0.6874 | 0.7144 | 0.0001 | 0.7239 | 0.8260 |
| 0.8984 | 0.9756 | 0.5833 | 0.5153 | 0.6119 | 0.0001 | 0.9845 | 0.3890 |
| 0.8991 | 0.7482 | 0.8085 | 0.6182 | 0.5584 | 0.0000 | 0.9063 | 0.8219 |
| 0.9845 | 0.3077 | 0.7060 | 0.3205 | 0.6976 | 0.0000 | 0.8651 | 0.3183 |
| 0.9273 | 0.8198 | 0.6673 | 0.2662 | 0.6263 | 0.0000 | 0.8552 | 0.9355 |
| 0.9408 | 0.8491 | 0.5793 | 0.0639 | 0.7034 | 0.0000 | 0.8048 | 0.6104 |
| 0.5067 | 0.9064 | 0.4457 | 0.5131 | 0.7748 | 0.0002 | 0.6702 | 0.2103 |
| 0.9497 | 0.7820 | 0.6386 | 0.2129 | 0.7174 | 0.0000 | 0.7313 | 0.9101 |
| 0.8641 | 0.5037 | 0.7951 | 0.3243 | 0.7143 | 0.0000 | 0.8313 | 0.1307 |
| 0.9226 | 0.5577 | 0.2943 | 0.7367 | 0.6620 | 0.0000 | 0.8858 | 0.4458 |
| 0.8184 | 0.7489 | 0.3294 | 0.1612 | 0.6157 | 0.0002 | 0.9570 | 0.9694 |
| 0.8376 | 0.8957 | 0.8602 | 0.3685 | 0.6976 | 0.0000 | 0.8293 | 0.7249 |
| 0.8568 | 0.6273 | 0.7714 | 0.3999 | 0.5505 | 0.0000 | 0.9528 | 0.9163 |
| 0.9257 | 0.8695 | 0.8274 | 0.2493 | 0.5280 | 0.0000 | 0.9722 | 0.4915 |
| 0.7676 | 0.7329 | 0.3333 | 0.7593 | 0.6514 | 0.0003 | 0.6667 | 0.2223 |
| 0.9568 | 0.6017 | 0.6134 | 0.8543 | 0.7344 | 0.0001 | 0.7741 | 0.2379 |
| 0.9352 | 0.7278 | 0.7917 | 0.5233 | 0.5331 | 0.0000 | 0.8746 | 0.4344 |
| 0.9355 | 0.6916 | 0.6623 | 0.6679 | 0.6615 | 0.0000 | 0.9036 | 0.0183 |
| 0.6310 | 0.8502 | 0.6220 | 0.3562 | 0.5929 | 0.0003 | 0.7791 | 0.3228 |
| 0.8967 | 0.7209 | 0.6902 | 0.5965 | 0.6511 | 0.0001 | 0.9359 | 0.6138 |
| 0.9649 | 0.3376 | 0.9251 | 0.3360 | 0.3971 | 0.0000 | 0.6356 | 0.3326 |
| 0.7165 | 0.6070 | 0.4904 | 0.5929 | 0.5743 | 0.0000 | 0.8057 | 0.6893 |
| 0.9065 | 0.8757 | 0.9212 | 0.3427 | 0.5647 | 0.0001 | 0.9682 | 0.0108 |
| 0.9419 | 0.7135 | 0.5563 | 0.7603 | 0.5819 | 0.0000 | 0.9792 | 0.4039 |
| 0.8797 | 0.8866 | 0.8290 | 0.8128 | 0.7100 | 0.0000 | 0.8553 | 0.6013 |
| 0.9253 | 0.7555 | 0.8126 | 0.7280 | 0.7534 | 0.0000 | 0.5947 | 0.7258 |
| 0.8570 | 0.6061 | 0.6241 | 0.1877 | 0.4910 | 0.0000 | 0.9774 | 0.9982 |
| 0.8948 | 0.7237 | 0.6609 | 0.7793 | 0.5435 | 0.0000 | 0.7200 | 0.6252 |
| 0.7520 | 0.8284 | 0.3779 | 0.8323 | 0.5898 | 0.0005 | 0.9139 | 0.1942 |
| 0.8961 | 0.8264 | 0.6969 | 0.4435 | 0.5193 | 0.0000 | 0.9700 | 0.0466 |
| 0.9650 | 0.6603 | 0.8297 | 0.8993 | 0.6370 | 0.0000 | 0.9632 | 0.8061 |
| 0.8140 | 0.9410 | 0.3537 | 0.3348 | 0.6339 | 0.0000 | 0.9765 | 0.1890 |
| 0.9211 | 0.8180 | 0.9026 | 0.1879 | 0.7011 | 0.0003 | 0.9342 | 0.1555 |
| 0.8864 | 0.8444 | 0.9207 | 0.5544 | 0.6588 | 0.0001 | 0.7320 | 0.4940 |
| 0.9748 | 0.6340 | 0.4682 | 0.5109 | 0.6159 | 0.0000 | 0.8498 | 0.8406 |
| 0.7481 | 0.9398 | 0.6305 | 0.4639 | 0.6586 | 0.0000 | 0.9426 | 0.2418 |
| 0.8464 | 0.7014 | 0.8176 | 0.8926 | 0.7167 | 0.0000 | 0.7849 | 0.3339 |
| 0.8505 | 0.5389 | 0.8874 | 0.7617 | 0.5216 | 0.0000 | 0.9971 | 0.6367 |
| 0.8261 | 0.5824 | 0.0849 | 0.8595 | 0.6180 | 0.0000 | 0.5909 | 0.0180 |
| 0.9315 | 0.8716 | 0.8889 | 0.4501 | 0.5721 | 0.0000 | 0.9955 | 0.8440 |
| 0.8405 | 0.3275 | 0.6956 | 0.6717 | 0.5399 | 0.0002 | 0.9769 | 0.4578 |
| 0.7579 | 0.7462 | 0.6093 | 0.1044 | 0.6643 | 0.0002 | 0.7885 | 0.3603 |
| 0.9783 | 0.8062 | 0.5647 | 0.7555 | 0.7506 | 0.0000 | 0.6389 | 0.1909 |
| 0.9188 | 0.9179 | 0.8325 | 0.0346 | 0.5705 | 0.0000 | 0.8328 | 0.3735 |
| 0.7256 | 0.9482 | 0.6941 | 0.4809 | 0.5595 | 0.0002 | 0.9554 | 0.0879 |
| 0.9599 | 0.8962 | 0.6641 | 0.5226 | 0.6366 | 0.0000 | 0.9090 | 0.7793 |
| 0.9264 | 0.9333 | 0.8553 | 0.4966 | 0.6182 | 0.0001 | 0.5884 | 0.8800 |
| 0.8316 | 0.8596 | 0.6672 | 0.5482 | 0.5138 | 0.0001 | 0.8758 | 0.1128 |
| 0.8811 | 0.8205 | 0.4934 | 0.7500 | 0.6666 | 0.0000 | 0.7979 | 0.8001 |
| 0.8315 | 0.8309 | 0.6367 | 0.4686 | 0.5233 | 0.0005 | 0.9732 | 0.5715 |
| 0.9157 | 0.5903 | 0.8866 | 0.7272 | 0.5167 | 0.0000 | 0.8593 | 0.0241 |
| 0.8939 | 0.8740 | 0.7923 | 0.8578 | 0.5689 | 0.0001 | 0.9065 | 0.5580 |
| 0.8996 | 0.6473 | 0.3029 | 0.4205 | 0.6720 | 0.0000 | 0.5176 | 0.3141 |
| 0.8797 | 0.8447 | 0.4068 | 0.5682 | 0.5626 | 0.0005 | 0.9830 | 0.4413 |
| 0.9563 | 0.7948 | 0.4523 | 0.5873 | 0.6635 | 0.0002 | 0.9628 | 0.0772 |
| 0.8796 | 0.8418 | 0.8911 | 0.3660 | 0.5529 | 0.0001 | 0.8734 | 0.8695 |
| 0.9199 | 0.8504 | 0.8659 | 0.7934 | 0.5900 | 0.0000 | 0.9483 | 0.7669 |
| 0.6757 | 0.9397 | 0.4749 | 0.7670 | 0.4964 | 0.0000 | 0.9187 | 0.0781 |
| 0.6984 | 0.8184 | 0.3807 | 0.7044 | 0.6069 | 0.0000 | 0.9174 | 0.7738 |
| 0.6854 | 0.7467 | 0.5722 | 0.6830 | 0.5659 | 0.0009 | 0.9760 | 0.1754 |
| 0.8061 | 0.5993 | 0.5461 | 0.5873 | 0.3813 | 0.0000 | 0.7793 | 0.7960 |
| 0.8244 | 0.4599 | 0.7466 | 0.9248 | 0.5250 | 0.0001 | 0.9599 | 0.1678 |
| 0.9692 | 0.7943 | 0.7972 | 0.5277 | 0.5865 | 0.0000 | 0.8463 | 0.9771 |
| 0.9110 | 0.7911 | 0.7236 | 0.6691 | 0.5038 | 0.0000 | 0.9767 | 0.5627 |
| 0.9038 | 0.7399 | 0.5783 | 0.9800 | 0.6489 | 0.0001 | 0.9600 | 0.4068 |
| 0.9407 | 0.7748 | 0.7027 | 0.8759 | 0.7645 | 0.0000 | 0.5809 | 0.4401 |
| 0.7041 | 0.6403 | 0.6907 | 0.3354 | 0.5401 | 0.0000 | 0.9993 | 0.0154 |
| 0.7274 | 0.8273 | 0.8318 | 0.1702 | 0.6512 | 0.0002 | 0.7980 | 0.5694 |
| 0.8079 | 0.9552 | 0.7781 | 0.3397 | 0.6883 | 0.0000 | 0.8111 | 0.0796 |
| 0.8871 | 0.8234 | 0.6466 | 0.5574 | 0.7356 | 0.0004 | 0.8568 | 0.4433 |
| 0.9026 | 0.7273 | 0.4961 | 0.2641 | 0.5074 | 0.0002 | 0.8180 | 0.8059 |
| 0.9146 | 0.6424 | 0.7511 | 0.3371 | 0.4454 | 0.0000 | 0.8507 | 0.3423 |
| 0.9788 | 0.7326 | 0.6124 | 0.6991 | 0.5420 | 0.0000 | 0.8028 | 0.5779 |
| 0.8751 | 0.9327 | 0.8492 | 0.5410 | 0.5588 | 0.0000 | 0.8957 | 0.9971 |
| 0.9709 | 0.6728 | 0.8741 | 0.7395 | 0.7587 | 0.0001 | 0.6279 | 0.7047 |
| 0.9250 | 0.3420 | 0.9191 | 0.5327 | 0.7307 | 0.0001 | 0.6933 | 0.3226 |
| 0.9426 | 0.6735 | 0.6626 | 0.6961 | 0.6313 | 0.0003 | 0.8207 | 0.3363 |
| 0.6254 | 0.8596 | 0.5991 | 0.3226 | 0.6036 | 0.0000 | 0.9906 | 0.1000 |
| 0.7412 | 0.8029 | 0.6652 | 0.6529 | 0.5803 | 0.0001 | 0.9072 | 0.3143 |
| 0.8055 | 0.7036 | 0.5488 | 0.2372 | 0.6988 | 0.0000 | 0.9383 | 0.4437 |
| 0.8727 | 0.8420 | 0.5956 | 0.8329 | 0.5079 | 0.0008 | 0.8081 | 0.4157 |
| 0.8418 | 0.6701 | 0.8348 | 0.1093 | 0.7787 | 0.0000 | 0.7896 | 0.9457 |
| 0.9130 | 0.6324 | 0.8930 | 0.4830 | 0.5885 | 0.0001 | 0.9906 | 0.4929 |
| 0.9227 | 0.7344 | 0.8294 | 0.6428 | 0.6438 | 0.0001 | 0.9725 | 0.7807 |
| 0.7838 | 0.5815 | 0.1003 | 0.3870 | 0.6456 | 0.0000 | 0.6877 | 0.7107 |
| 0.8913 | 0.7154 | 0.8755 | 0.0469 | 0.7830 | 0.0000 | 0.1891 | 0.2918 |
| 0.9074 | 0.5843 | 0.8031 | 0.9347 | 0.7343 | 0.0000 | 0.8610 | 0.7805 |
| 0.9246 | 0.8589 | 0.5513 | 0.2134 | 0.6590 | 0.0000 | 0.8822 | 0.6235 |
| 0.8515 | 0.7611 | 0.5433 | 0.5983 | 0.5594 | 0.0000 | 0.8174 | 0.1758 |
| 0.8018 | 0.8284 | 0.4908 | 0.6346 | 0.7348 | 0.0000 | 0.8457 | 0.5693 |
| 0.5553 | 0.7644 | 0.3421 | 0.8564 | 0.4778 | 0.0000 | 0.9490 | 0.4400 |
| 0.9638 | 0.6715 | 0.5539 | 0.0040 | 0.3597 | 0.0009 | 0.9177 | 0.1076 |
| 0.8227 | 0.7183 | 0.3700 | 0.8331 | 0.6175 | 0.0000 | 0.8160 | 0.2104 |
| 0.8825 | 0.9423 | 0.7976 | 0.0811 | 0.4518 | 0.0001 | 0.9183 | 0.7242 |
| 0.6776 | 0.8638 | 0.7140 | 0.6558 | 0.5961 | 0.0000 | 0.8653 | 0.4989 |
| 0.9366 | 0.8894 | 0.8449 | 0.5546 | 0.6727 | 0.0000 | 0.9151 | 0.8341 |
| 0.6884 | 0.7858 | 0.5066 | 0.3641 | 0.6104 | 0.0006 | 0.9979 | 0.2015 |
| 0.9118 | 0.6912 | 0.8171 | 0.1605 | 0.6524 | 0.0001 | 0.9753 | 0.5807 |
| 0.7984 | 0.6920 | 0.8605 | 0.2831 | 0.5223 | 0.0000 | 0.9284 | 0.8706 |
| 0.9853 | 0.6921 | 0.2908 | 0.2587 | 0.6252 | 0.0000 | 0.6311 | 0.0480 |
| 0.9208 | 0.5721 | 0.5433 | 0.2869 | 0.5132 | 0.0000 | 0.8762 | 0.9765 |
| 0.7717 | 0.8495 | 0.7085 | 0.7425 | 0.6156 | 0.0000 | 0.9805 | 0.8348 |
| 0.9106 | 0.6234 | 0.9329 | 0.4635 | 0.7205 | 0.0000 | 0.8332 | 0.6800 |
| 0.9784 | 0.9105 | 0.7218 | 0.8199 | 0.5212 | 0.0000 | 0.9816 | 0.1629 |
| 0.7810 | 0.7513 | 0.7675 | 0.8345 | 0.5799 | 0.0000 | 0.9873 | 0.4096 |
| 0.8600 | 0.7790 | 0.2713 | 0.7506 | 0.5620 | 0.0000 | 0.3975 | 0.9990 |
| 0.8187 | 0.8181 | 0.8362 | 0.7775 | 0.5804 | 0.0000 | 0.8865 | 0.4697 |
| 0.8992 | 0.7663 | 0.6643 | 0.7660 | 0.4275 | 0.0000 | 0.8069 | 0.0240 |
| 0.5946 | 0.6773 | 0.1491 | 0.7667 | 0.5354 | 0.0000 | 0.9931 | 0.8858 |
| 0.9387 | 0.8886 | 0.8144 | 0.1518 | 0.5571 | 0.0001 | 0.9287 | 0.3807 |
| 0.9554 | 0.8190 | 0.4800 | 0.2893 | 0.4428 | 0.0000 | 0.9834 | 0.1441 |
| 0.9035 | 0.5085 | 0.6147 | 0.7003 | 0.7582 | 0.0000 | 0.8010 | 0.7420 |
| 0.8615 | 0.8709 | 0.4453 | 0.5726 | 0.6181 | 0.0000 | 0.8033 | 0.3200 |
| 0.7847 | 0.8347 | 0.7286 | 0.4975 | 0.6012 | 0.0000 | 0.9404 | 0.6620 |
| 0.7964 | 0.7174 | 0.1102 | 0.6733 | 0.6107 | 0.0000 | 0.9882 | 0.2082 |
| 0.9728 | 0.8177 | 0.8267 | 0.5714 | 0.6559 | 0.0002 | 0.9583 | 0.7412 |
| 0.4279 | 0.8401 | 0.4629 | 0.1035 | 0.4903 | 0.0001 | 0.7356 | 0.6052 |
| 0.8580 | 0.6386 | 0.8321 | 0.6505 | 0.4369 | 0.0000 | 0.8945 | 0.5873 |
| 0.9011 | 0.5253 | 0.7119 | 0.3660 | 0.7121 | 0.0000 | 0.8478 | 0.3937 |
| 0.9164 | 0.8681 | 0.6796 | 0.2437 | 0.6364 | 0.0000 | 0.7813 | 0.5822 |
| 0.9512 | 0.7872 | 0.9051 | 0.5246 | 0.5297 | 0.0001 | 0.9999 | 0.1227 |
| 0.9549 | 0.4548 | 0.1899 | 0.7806 | 0.7683 | 0.0003 | 0.2535 | 0.7228 |
| 0.9703 | 0.7399 | 0.4941 | 0.6435 | 0.6498 | 0.0010 | 0.8241 | 0.4164 |
| 0.6195 | 0.9045 | 0.6744 | 0.2606 | 0.5971 | 0.0001 | 0.9908 | 0.4249 |
| 0.9603 | 0.4299 | 0.7246 | 0.4418 | 0.5554 | 0.0001 | 0.9223 | 0.3902 |
| 0.5718 | 0.6204 | 0.5571 | 0.7984 | 0.5490 | 0.0000 | 0.9160 | 0.0699 |
| 0.7186 | 0.7443 | 0.4456 | 0.2862 | 0.5962 | 0.0001 | 0.9255 | 0.1765 |
| 0.7899 | 0.8747 | 0.2320 | 0.2934 | 0.5885 | 0.0000 | 0.9818 | 0.9682 |
| 0.8936 | 0.5313 | 0.3492 | 0.1740 | 0.6306 | 0.0000 | 0.8642 | 0.1498 |
| 0.8800 | 0.9335 | 0.8462 | 0.4695 | 0.5599 | 0.0002 | 0.8689 | 0.3158 |
| 0.9249 | 0.7601 | 0.6880 | 0.0694 | 0.6076 | 0.0000 | 0.9492 | 0.3096 |
| 0.9364 | 0.9549 | 0.8269 | 0.6417 | 0.7195 | 0.0000 | 0.6568 | 0.9887 |
| 0.7267 | 0.7493 | 0.6239 | 0.1427 | 0.6808 | 0.0000 | 0.9753 | 0.9820 |
| 0.8179 | 0.7567 | 0.6534 | 0.5403 | 0.5982 | 0.0000 | 0.9721 | 0.4330 |
| 0.7792 | 0.8198 | 0.8348 | 0.3562 | 0.6057 | 0.0000 | 0.9651 | 0.9044 |
| 0.8628 | 0.8621 | 0.7779 | 0.1810 | 0.6377 | 0.0000 | 0.9759 | 0.8835 |
| 0.7993 | 0.9007 | 0.5915 | 0.6654 | 0.6106 | 0.0000 | 0.9898 | 0.3542 |
| 0.8613 | 0.9058 | 0.9110 | 0.1116 | 0.5000 | 0.0000 | 0.9638 | 0.6659 |
| 0.7982 | 0.8347 | 0.8834 | 0.4879 | 0.6477 | 0.0000 | 0.8793 | 0.5449 |
| 0.8109 | 0.6122 | 0.7193 | 0.9656 | 0.6481 | 0.0000 | 0.9654 | 0.1257 |
| 0.9182 | 0.8047 | 0.6586 | 0.1408 | 0.6439 | 0.0000 | 0.7063 | 0.6776 |
| 0.7961 | 0.8760 | 0.3765 | 0.5330 | 0.5339 | 0.0000 | 0.9830 | 0.9720 |
| 0.8885 | 0.6685 | 0.6298 | 0.4869 | 0.5966 | 0.0001 | 0.9947 | 0.7716 |
| 0.7966 | 0.8140 | 0.6298 | 0.1299 | 0.6227 | 0.0000 | 0.9514 | 0.9449 |
| 0.8916 | 0.9302 | 0.4642 | 0.2785 | 0.7088 | 0.0005 | 0.8381 | 0.8466 |
| 0.9150 | 0.8940 | 0.6833 | 0.5760 | 0.6626 | 0.0001 | 0.8651 | 0.6991 |
| 0.9667 | 0.5423 | 0.4158 | 0.3536 | 0.7870 | 0.0000 | 0.6281 | 0.4522 |
| 0.9445 | 0.6627 | 0.9157 | 0.5750 | 0.7101 | 0.0000 | 0.8159 | 0.9900 |
| 0.6854 | 0.7517 | 0.3394 | 0.5792 | 0.5470 | 0.0000 | 0.8831 | 0.8135 |
| 0.7637 | 0.6932 | 0.5937 | 0.2385 | 0.6245 | 0.0000 | 0.7614 | 0.6523 |
| 0.9719 | 0.6126 | 0.4774 | 0.6116 | 0.6228 | 0.0000 | 0.8638 | 0.9937 |
| 0.8978 | 0.6114 | 0.6367 | 0.4403 | 0.6142 | 0.0000 | 0.9028 | 0.5676 |
| 0.9232 | 0.5316 | 0.7988 | 0.9372 | 0.8447 | 0.0000 | 0.2847 | 0.8447 |
| 0.9067 | 0.6137 | 0.9054 | 0.8033 | 0.6130 | 0.0008 | 0.9698 | 0.3404 |
| 0.9502 | 0.9324 | 0.6573 | 0.7400 | 0.6560 | 0.0003 | 0.8760 | 0.5575 |
| 0.9054 | 0.4577 | 0.7581 | 0.3817 | 0.5077 | 0.0000 | 0.7707 | 0.8000 |
| 0.9437 | 0.6178 | 0.9063 | 0.1743 | 0.5792 | 0.0001 | 0.8738 | 0.5363 |
| 0.8382 | 0.7657 | 0.1563 | 0.7613 | 0.7357 | 0.0000 | 0.5152 | 0.6569 |
| 0.9665 | 0.9772 | 0.6654 | 0.3296 | 0.7027 | 0.0001 | 0.8402 | 0.7632 |
| 0.8760 | 0.6405 | 0.8626 | 0.9135 | 0.5878 | 0.0000 | 0.8268 | 0.9894 |
| 0.9442 | 0.4949 | 0.8226 | 0.0352 | 0.6954 | 0.0000 | 0.6198 | 0.9248 |
| 0.8942 | 0.9125 | 0.7574 | 0.5114 | 0.5360 | 0.0005 | 0.9875 | 0.9748 |
| 0.9683 | 0.8997 | 0.5938 | 0.7391 | 0.5766 | 0.0000 | 0.9560 | 0.5499 |
| 0.9111 | 0.8925 | 0.7238 | 0.6316 | 0.5604 | 0.0000 | 0.8447 | 0.6604 |
| 0.7981 | 0.8187 | 0.3224 | 0.6360 | 0.7287 | 0.0000 | 0.9251 | 0.6836 |
| 0.9103 | 0.3460 | 0.6437 | 0.1332 | 0.7443 | 0.0001 | 0.7443 | 0.2825 |
| 0.8446 | 0.5172 | 0.5656 | 0.7787 | 0.6434 | 0.0001 | 0.9190 | 0.0740 |
| 0.7879 | 0.6559 | 0.7035 | 0.6676 | 0.5248 | 0.0012 | 0.9720 | 0.0163 |
| 0.9052 | 0.5499 | 0.7401 | 0.1044 | 0.6890 | 0.0000 | 0.7149 | 0.6582 |
| 0.8572 | 0.7133 | 0.4742 | 0.4272 | 0.5908 | 0.0003 | 0.9076 | 0.2919 |
| 0.8995 | 0.3747 | 0.7673 | 0.6065 | 0.5305 | 0.0000 | 0.9530 | 0.7589 |
| 0.7921 | 0.8591 | 0.2503 | 0.4773 | 0.6326 | 0.0000 | 0.5836 | 0.9647 |
| 0.8887 | 0.6260 | 0.7931 | 0.7350 | 0.6875 | 0.0000 | 0.9360 | 0.3311 |
| 0.8919 | 0.8068 | 0.8600 | 0.7577 | 0.6922 | 0.0001 | 0.8671 | 0.4295 |
| 0.8761 | 0.4324 | 0.6243 | 0.7790 | 0.7509 | 0.0000 | 0.4216 | 0.7451 |
| 0.8968 | 0.9250 | 0.1071 | 0.7947 | 0.8429 | 0.0000 | 0.2275 | 0.3481 |
| 0.9096 | 0.2940 | 0.7766 | 0.9673 | 0.7050 | 0.0000 | 0.8724 | 0.4502 |
| 0.6807 | 0.7559 | 0.1396 | 0.8057 | 0.7409 | 0.0001 | 0.4544 | 0.1427 |
| 0.9152 | 0.6365 | 0.7562 | 0.1814 | 0.5643 | 0.0000 | 0.9065 | 0.6587 |
| 0.9197 | 0.8811 | 0.6345 | 0.6748 | 0.7064 | 0.0000 | 0.8659 | 0.8702 |
| 0.9765 | 0.6545 | 0.7542 | 0.8859 | 0.6589 | 0.0000 | 0.8112 | 0.9746 |
| 0.9388 | 0.9867 | 0.9364 | 0.2342 | 0.7701 | 0.0000 | 0.8783 | 0.4862 |
| 0.6774 | 0.7828 | 0.5588 | 0.8019 | 0.5672 | 0.0001 | 0.8632 | 0.4258 |
| 0.7976 | 0.8292 | 0.7267 | 0.0149 | 0.5333 | 0.0000 | 0.9565 | 0.4772 |
| 0.6052 | 0.8332 | 0.3808 | 0.6805 | 0.5684 | 0.0000 | 0.9721 | 0.7022 |
| 0.9231 | 0.6461 | 0.8807 | 0.4116 | 0.6754 | 0.0000 | 0.9397 | 0.9453 |
| 0.8693 | 0.9223 | 0.7930 | 0.3624 | 0.7510 | 0.0000 | 0.3892 | 0.4280 |
| 0.9475 | 0.5520 | 0.6060 | 0.2246 | 0.7243 | 0.0000 | 0.7470 | 0.9710 |
| 0.8501 | 0.4378 | 0.5937 | 0.4780 | 0.6206 | 0.0000 | 0.9645 | 0.1608 |
| 0.8351 | 0.9265 | 0.7757 | 0.3280 | 0.5056 | 0.0000 | 0.9969 | 0.8628 |
| 0.8897 | 0.8044 | 0.8567 | 0.3786 | 0.7358 | 0.0000 | 0.8639 | 0.9623 |
| 0.8361 | 0.4250 | 0.7213 | 0.2971 | 0.5397 | 0.0001 | 0.9738 | 0.0446 |
| 0.7895 | 0.7104 | 0.5824 | 0.3479 | 0.5342 | 0.0000 | 0.8434 | 0.5155 |
| 0.8214 | 0.9879 | 0.6134 | 0.5978 | 0.6698 | 0.0001 | 0.8367 | 0.7203 |
| 0.8928 | 0.6191 | 0.7353 | 0.9286 | 0.4939 | 0.0000 | 0.8313 | 0.8910 |
| 0.9091 | 0.8389 | 0.6997 | 0.1537 | 0.6871 | 0.0000 | 0.7266 | 0.7523 |
| 0.9205 | 0.6643 | 0.8289 | 0.8554 | 0.7827 | 0.0000 | 0.5410 | 0.6974 |
| 0.7701 | 0.8408 | 0.5692 | 0.7162 | 0.6289 | 0.0000 | 0.9746 | 0.9196 |
| 0.6583 | 0.5752 | 0.6298 | 0.3543 | 0.6514 | 0.0000 | 0.7358 | 0.7904 |
| 0.9046 | 0.9399 | 0.6006 | 0.4357 | 0.7782 | 0.0000 | 0.8943 | 0.6537 |
| 0.7603 | 0.5022 | 0.3719 | 0.3199 | 0.7066 | 0.0000 | 0.7994 | 0.8272 |
| 0.8198 | 0.5961 | 0.5904 | 0.8281 | 0.4626 | 0.0000 | 0.9693 | 0.9486 |
| 0.8698 | 0.7707 | 0.7424 | 0.4192 | 0.6293 | 0.0000 | 0.8170 | 0.4581 |
| 0.3600 | 0.7621 | 0.2056 | 0.9267 | 0.5270 | 0.0000 | 0.9994 | 0.3343 |
| 0.9493 | 0.7493 | 0.6184 | 0.8754 | 0.5051 | 0.0001 | 0.8047 | 0.0701 |
| 0.8133 | 0.8980 | 0.6321 | 0.3181 | 0.7436 | 0.0000 | 0.8526 | 0.9774 |
| 0.8674 | 0.6621 | 0.6922 | 0.6581 | 0.6496 | 0.0003 | 0.9670 | 0.7357 |
| 0.9218 | 0.8516 | 0.7465 | 0.3698 | 0.6993 | 0.0000 | 0.8918 | 0.3595 |
| 0.7987 | 0.7143 | 0.7962 | 0.6442 | 0.6208 | 0.0001 | 0.7877 | 0.0520 |
| 0.9602 | 0.9556 | 0.8124 | 0.6345 | 0.6600 | 0.0003 | 0.9422 | 0.0178 |
| 0.8963 | 0.7918 | 0.6214 | 0.5306 | 0.3973 | 0.0000 | 0.9921 | 0.8879 |
| 0.6816 | 0.5997 | 0.6986 | 0.5923 | 0.5097 | 0.0000 | 0.9340 | 0.5608 |
| 0.8808 | 0.4588 | 0.7435 | 0.3888 | 0.7169 | 0.0000 | 0.9168 | 0.5985 |
| 0.9641 | 0.8088 | 0.7729 | 0.1993 | 0.5318 | 0.0000 | 0.9830 | 0.7179 |
| 0.7815 | 0.8445 | 0.5484 | 0.4170 | 0.4845 | 0.0018 | 0.9571 | 0.1839 |
| 0.8902 | 0.5398 | 0.6172 | 0.4282 | 0.5878 | 0.0002 | 0.9749 | 0.1105 |
| 0.9116 | 0.8114 | 0.4044 | 0.6529 | 0.6779 | 0.0004 | 0.7002 | 0.1531 |
| 0.9366 | 0.8295 | 0.4414 | 0.5696 | 0.7785 | 0.0000 | 0.7929 | 0.6003 |
| 0.9493 | 0.7944 | 0.8538 | 0.6321 | 0.7307 | 0.0002 | 0.7909 | 0.1291 |
| 0.8887 | 0.8463 | 0.4764 | 0.3400 | 0.5197 | 0.0003 | 0.9060 | 0.9581 |
| 0.8374 | 0.7963 | 0.4453 | 0.6901 | 0.6838 | 0.0000 | 0.7256 | 0.7352 |
| 0.8755 | 0.8644 | 0.7790 | 0.6892 | 0.7068 | 0.0001 | 0.8632 | 0.9111 |
| 0.6281 | 0.9058 | 0.1835 | 0.6820 | 0.6629 | 0.0000 | 0.8789 | 0.1158 |
| 0.9524 | 0.8468 | 0.3983 | 0.8069 | 0.6438 | 0.0000 | 0.4023 | 0.6004 |
| 0.9069 | 0.5446 | 0.5626 | 0.9790 | 0.7201 | 0.0000 | 0.7652 | 0.9729 |
| 0.8900 | 0.7895 | 0.7742 | 0.8029 | 0.7388 | 0.0001 | 0.7563 | 0.4627 |
| 0.9584 | 0.8160 | 0.8111 | 0.4469 | 0.6290 | 0.0001 | 0.8831 | 0.3743 |
| 0.9131 | 0.5272 | 0.6363 | 0.8420 | 0.5193 | 0.0001 | 0.9932 | 0.7225 |
| 0.7863 | 0.8329 | 0.5726 | 0.0680 | 0.8120 | 0.0001 | 0.1599 | 0.1734 |
| 0.8820 | 0.9324 | 0.7068 | 0.7731 | 0.5605 | 0.0001 | 0.7303 | 0.6290 |
| 0.8359 | 0.3932 | 0.5764 | 0.3894 | 0.4748 | 0.0001 | 0.9432 | 0.4061 |
| 0.9658 | 0.6492 | 0.5512 | 0.6548 | 0.6326 | 0.0000 | 0.8974 | 0.8112 |
| 0.8553 | 0.7792 | 0.6559 | 0.5046 | 0.7053 | 0.0003 | 0.7058 | 0.0829 |
| 0.8051 | 0.8311 | 0.4977 | 0.2118 | 0.5952 | 0.0002 | 0.9601 | 0.5788 |
| 0.9270 | 0.8762 | 0.7263 | 0.3780 | 0.6615 | 0.0000 | 0.9637 | 0.2075 |
| 0.9459 | 0.5951 | 0.6867 | 0.3860 | 0.5218 | 0.0000 | 0.8378 | 0.9730 |
| 0.7655 | 0.7752 | 0.6894 | 0.2040 | 0.6807 | 0.0000 | 0.7505 | 0.3286 |
| 0.8682 | 0.8563 | 0.8983 | 0.7384 | 0.5816 | 0.0000 | 0.9860 | 0.9721 |
| 0.9484 | 0.8011 | 0.7394 | 0.7871 | 0.4246 | 0.0001 | 0.9820 | 0.3215 |
| 0.6986 | 0.7717 | 0.4926 | 0.2684 | 0.6703 | 0.0001 | 0.9057 | 0.7608 |
| 0.9420 | 0.5258 | 0.7695 | 0.9099 | 0.4882 | 0.0000 | 0.8133 | 0.8260 |
| 0.8855 | 0.8077 | 0.6122 | 0.3130 | 0.4989 | 0.0000 | 0.9986 | 0.8187 |
| 0.8422 | 0.7407 | 0.4291 | 0.1988 | 0.6417 | 0.0000 | 0.8785 | 0.9334 |
| 0.7655 | 0.5505 | 0.7807 | 0.3849 | 0.6439 | 0.0000 | 0.9710 | 0.9981 |
| 0.9404 | 0.3690 | 0.7891 | 0.5555 | 0.5665 | 0.0000 | 0.8482 | 0.9970 |
| 0.6000 | 0.7984 | 0.3827 | 0.9130 | 0.6529 | 0.0000 | 0.8964 | 0.5541 |
| 0.9259 | 0.8062 | 0.3689 | 0.5570 | 0.6588 | 0.0002 | 0.9868 | 0.8666 |
| 0.6626 | 0.9747 | 0.3079 | 0.1784 | 0.6438 | 0.0000 | 0.8301 | 0.3266 |
| 0.7747 | 0.5491 | 0.7020 | 0.4761 | 0.5573 | 0.0001 | 0.8130 | 0.8746 |
| 0.9591 | 0.7777 | 0.7928 | 0.2087 | 0.5953 | 0.0002 | 0.4154 | 0.9665 |
| 0.9823 | 0.6352 | 0.6764 | 0.4782 | 0.6440 | 0.0001 | 0.8453 | 0.3249 |
| 0.8390 | 0.6407 | 0.4781 | 0.6333 | 0.5362 | 0.0002 | 0.9736 | 0.2544 |
| 0.8661 | 0.6701 | 0.8915 | 0.7868 | 0.4836 | 0.0000 | 0.8508 | 0.4172 |
| 0.9077 | 0.7470 | 0.5857 | 0.0441 | 0.5433 | 0.0000 | 0.8739 | 0.5374 |
| 0.8644 | 0.6502 | 0.1785 | 0.7812 | 0.7084 | 0.0000 | 0.7854 | 0.2742 |
| 0.9574 | 0.9297 | 0.7259 | 0.5338 | 0.7582 | 0.0001 | 0.6365 | 0.0072 |
| 0.8384 | 0.8832 | 0.4027 | 0.0424 | 0.6436 | 0.0002 | 0.8359 | 0.6685 |
| 0.8598 | 0.8129 | 0.8404 | 0.7320 | 0.5418 | 0.0000 | 0.7862 | 0.6162 |
| 0.9867 | 0.7887 | 0.7447 | 0.6200 | 0.6244 | 0.0000 | 0.6977 | 0.1644 |
| 0.8425 | 0.7024 | 0.6184 | 0.2705 | 0.6377 | 0.0001 | 0.8624 | 0.0767 |
| 0.9333 | 0.8187 | 0.6483 | 0.4982 | 0.4757 | 0.0001 | 0.9245 | 0.0195 |
| 0.9515 | 0.8864 | 0.4489 | 0.1834 | 0.6689 | 0.0001 | 0.9539 | 0.1875 |
| 0.9320 | 0.8795 | 0.6854 | 0.4639 | 0.7414 | 0.0000 | 0.3032 | 0.3564 |
| 0.8610 | 0.7328 | 0.6743 | 0.5473 | 0.6331 | 0.0000 | 0.5785 | 0.6895 |
| 0.8031 | 0.8834 | 0.6121 | 0.1626 | 0.5132 | 0.0000 | 0.9969 | 0.4733 |
| 0.9273 | 0.8583 | 0.8724 | 0.1911 | 0.4671 | 0.0001 | 0.8910 | 0.7910 |
| 0.9296 | 0.6611 | 0.5560 | 0.6594 | 0.6686 | 0.0001 | 0.7854 | 0.7589 |
| 0.8353 | 0.7015 | 0.6988 | 0.9453 | 0.5784 | 0.0000 | 0.9486 | 0.3402 |
| 0.9566 | 0.9054 | 0.5453 | 0.4460 | 0.6495 | 0.0000 | 0.9616 | 0.6186 |
| 0.9456 | 0.9172 | 0.6692 | 0.5072 | 0.6228 | 0.0001 | 0.8578 | 0.3896 |
| 0.8483 | 0.7614 | 0.2132 | 0.6971 | 0.7129 | 0.0002 | 0.7066 | 0.6023 |
| 0.7778 | 0.7270 | 0.4321 | 0.5129 | 0.5915 | 0.0000 | 0.9915 | 0.0349 |
| 0.8795 | 0.7861 | 0.5978 | 0.6807 | 0.6358 | 0.0001 | 0.9637 | 0.8594 |
| 0.8607 | 0.7138 | 0.5925 | 0.3767 | 0.5725 | 0.0000 | 0.9003 | 0.7220 |
| 0.8577 | 0.4875 | 0.8410 | 0.4749 | 0.6434 | 0.0000 | 0.5268 | 0.9820 |
| 0.9296 | 0.5183 | 0.8047 | 0.1823 | 0.7358 | 0.0008 | 0.7629 | 0.0427 |
| 0.9088 | 0.4647 | 0.8125 | 0.2943 | 0.6335 | 0.0000 | 0.4691 | 0.9290 |
| 0.9566 | 0.6878 | 0.8016 | 0.2813 | 0.6794 | 0.0000 | 0.9796 | 0.0657 |
| 0.9529 | 0.7168 | 0.3894 | 0.7270 | 0.6444 | 0.0000 | 0.6522 | 0.8300 |
| 0.9544 | 0.8370 | 0.6144 | 0.4827 | 0.7359 | 0.0000 | 0.8832 | 0.5154 |
| 0.8884 | 0.5618 | 0.2776 | 0.5173 | 0.7732 | 0.0000 | 0.4434 | 0.2613 |
| 0.7536 | 0.7391 | 0.3289 | 0.6444 | 0.5877 | 0.0000 | 0.9947 | 0.8715 |
| 0.7870 | 0.9022 | 0.7938 | 0.7471 | 0.5339 | 0.0000 | 0.7875 | 0.7648 |
| 0.8396 | 0.7896 | 0.8298 | 0.7621 | 0.5088 | 0.0006 | 0.9862 | 0.1490 |
| 0.8835 | 0.7017 | 0.6814 | 0.9837 | 0.7428 | 0.0004 | 0.8976 | 0.6827 |
| 0.8976 | 0.7709 | 0.8431 | 0.4380 | 0.6405 | 0.0000 | 0.9791 | 0.8488 |
| 0.9243 | 0.6003 | 0.7618 | 0.7296 | 0.5594 | 0.0000 | 0.8728 | 0.8150 |
| 0.9371 | 0.7297 | 0.4649 | 0.6705 | 0.7300 | 0.0001 | 0.5847 | 0.3424 |
| 0.9172 | 0.6349 | 0.8114 | 0.7824 | 0.5363 | 0.0000 | 0.9071 | 0.5942 |
| 0.9796 | 0.9099 | 0.5897 | 0.3707 | 0.7439 | 0.0001 | 0.7550 | 0.6592 |
| 0.7979 | 0.8970 | 0.4252 | 0.6925 | 0.4680 | 0.0000 | 0.7066 | 0.4730 |
| 0.9310 | 0.6934 | 0.5033 | 0.5322 | 0.4882 | 0.0000 | 0.7886 | 0.7945 |
| 0.9229 | 0.8137 | 0.8191 | 0.8334 | 0.5573 | 0.0001 | 0.7690 | 0.9540 |
| 0.9845 | 0.6760 | 0.5064 | 0.2138 | 0.7210 | 0.0000 | 0.7766 | 0.5599 |
| 0.8619 | 0.9429 | 0.5421 | 0.1834 | 0.5159 | 0.0004 | 0.9979 | 0.1016 |
| 0.9937 | 0.6724 | 0.1760 | 0.5907 | 0.7448 | 0.0003 | 0.3196 | 0.6085 |
| 0.9582 | 0.5019 | 0.8879 | 0.3259 | 0.5919 | 0.0000 | 0.9381 | 0.6350 |
| 0.9584 | 0.9415 | 0.9033 | 0.4295 | 0.4405 | 0.0003 | 0.9752 | 0.7599 |
| 0.9496 | 0.6106 | 0.6512 | 0.6738 | 0.6565 | 0.0000 | 0.9726 | 0.3309 |
| 0.8787 | 0.6919 | 0.3594 | 0.1683 | 0.5406 | 0.0001 | 0.9339 | 0.6447 |
| 0.8123 | 0.9837 | 0.6951 | 0.1363 | 0.6542 | 0.0000 | 0.7977 | 0.9030 |
| 0.5114 | 0.7216 | 0.4806 | 0.3782 | 0.6694 | 0.0000 | 0.9182 | 0.4823 |
| 0.8628 | 0.7440 | 0.7862 | 0.8497 | 0.5628 | 0.0001 | 0.8713 | 0.6892 |
| 0.9437 | 0.7188 | 0.7148 | 0.5992 | 0.6062 | 0.0000 | 0.7894 | 0.9363 |
| 0.8493 | 0.8851 | 0.5808 | 0.5172 | 0.6034 | 0.0003 | 0.9519 | 0.5364 |
| 0.6115 | 0.8915 | 0.4900 | 0.2888 | 0.5127 | 0.0000 | 0.9978 | 0.3809 |
| 0.9098 | 0.9470 | 0.7428 | 0.4215 | 0.6444 | 0.0000 | 0.6906 | 0.0060 |
| 0.9575 | 0.9643 | 0.6431 | 0.2007 | 0.5801 | 0.0000 | 0.8313 | 0.9778 |
| 0.7975 | 0.5409 | 0.8094 | 0.7263 | 0.7089 | 0.0000 | 0.9203 | 0.7037 |
| 0.9176 | 0.4555 | 0.6796 | 0.4371 | 0.6921 | 0.0000 | 0.7688 | 0.5143 |
| 0.8755 | 0.6139 | 0.5570 | 0.2895 | 0.6865 | 0.0000 | 0.7443 | 0.9685 |
| 0.8909 | 0.7549 | 0.7366 | 0.7679 | 0.6459 | 0.0000 | 0.9157 | 0.9973 |
| 0.9246 | 0.8668 | 0.7538 | 0.6139 | 0.7056 | 0.0000 | 0.8690 | 0.2458 |
| 0.9479 | 0.8403 | 0.8677 | 0.8965 | 0.7057 | 0.0000 | 0.6042 | 0.9965 |
| 0.9040 | 0.7571 | 0.6637 | 0.4005 | 0.5814 | 0.0001 | 0.6816 | 0.7546 |
| 0.9732 | 0.8357 | 0.1170 | 0.8291 | 0.6551 | 0.0000 | 0.9233 | 0.8737 |
| 0.9864 | 0.8916 | 0.5270 | 0.2640 | 0.5959 | 0.0000 | 0.7286 | 0.3881 |
| 0.9355 | 0.7289 | 0.5704 | 0.8080 | 0.7075 | 0.0000 | 0.8077 | 0.5169 |
| 0.9716 | 0.7762 | 0.7763 | 0.3576 | 0.6945 | 0.0001 | 0.7416 | 0.8107 |
| 0.9113 | 0.4589 | 0.4641 | 0.8254 | 0.8256 | 0.0000 | 0.5198 | 0.2824 |
| 0.8521 | 0.8824 | 0.6215 | 0.3647 | 0.5904 | 0.0000 | 0.8876 | 0.8120 |
| 0.8859 | 0.8119 | 0.2248 | 0.5825 | 0.6147 | 0.0000 | 0.9505 | 0.5973 |
| 0.8732 | 0.6845 | 0.4902 | 0.6853 | 0.6073 | 0.0000 | 0.9615 | 0.7657 |
| 0.8133 | 0.8965 | 0.4798 | 0.7754 | 0.5403 | 0.0001 | 0.9481 | 0.3409 |
| 0.8374 | 0.8679 | 0.7752 | 0.6796 | 0.5447 | 0.0000 | 0.7799 | 0.7650 |
| 0.9339 | 0.7514 | 0.7199 | 0.4464 | 0.6193 | 0.0000 | 0.9167 | 0.8216 |
| 0.9023 | 0.9342 | 0.8354 | 0.0846 | 0.4995 | 0.0001 | 1.0000 | 0.5601 |
| 0.9535 | 0.7329 | 0.4394 | 0.1220 | 0.6851 | 0.0002 | 0.9830 | 0.3955 |
| 0.9245 | 0.7181 | 0.8318 | 0.9360 | 0.5623 | 0.0000 | 0.8434 | 0.9458 |
| 0.8593 | 0.5987 | 0.8643 | 0.6970 | 0.7134 | 0.0000 | 0.6962 | 0.3599 |
| 0.9696 | 0.7195 | 0.7186 | 0.8421 | 0.4644 | 0.0001 | 0.8761 | 0.5616 |
| 0.9419 | 0.6489 | 0.7799 | 0.5058 | 0.7884 | 0.0001 | 0.2504 | 0.6870 |
| 0.8325 | 0.6898 | 0.8285 | 0.3579 | 0.7146 | 0.0001 | 0.6036 | 0.6878 |
| 0.7285 | 0.5247 | 0.4468 | 0.2989 | 0.5827 | 0.0003 | 0.9518 | 0.4747 |
| 0.7687 | 0.8309 | 0.4829 | 0.9075 | 0.6744 | 0.0000 | 0.8613 | 0.4640 |
| 0.8344 | 0.8294 | 0.6970 | 0.7812 | 0.4745 | 0.0000 | 0.9993 | 0.8000 |
| 0.8428 | 0.7978 | 0.5732 | 0.4995 | 0.6643 | 0.0010 | 0.7928 | 0.3513 |
| 0.9800 | 0.6598 | 0.6801 | 0.0720 | 0.7146 | 0.0002 | 0.6786 | 0.8652 |
| 0.8904 | 0.7004 | 0.6127 | 0.2811 | 0.6498 | 0.0001 | 0.7865 | 0.0453 |
| 0.9483 | 0.6010 | 0.5610 | 0.8518 | 0.5525 | 0.0000 | 0.8845 | 0.5594 |
| 0.9171 | 0.7596 | 0.5616 | 0.4224 | 0.5829 | 0.0000 | 0.8736 | 0.0741 |
| 0.9115 | 0.3569 | 0.6277 | 0.6330 | 0.6895 | 0.0000 | 0.8719 | 0.8076 |
| 0.8218 | 0.6180 | 0.7292 | 0.6647 | 0.5007 | 0.0021 | 0.9736 | 0.1752 |
| 0.7710 | 0.7420 | 0.5667 | 0.7188 | 0.7106 | 0.0000 | 0.7200 | 0.1083 |
| 0.8094 | 0.7982 | 0.4661 | 0.6292 | 0.6204 | 0.0001 | 0.9796 | 0.2381 |
| 0.8116 | 0.6701 | 0.5906 | 0.8226 | 0.6066 | 0.0000 | 0.9578 | 0.2396 |
| 0.8182 | 0.7162 | 0.7230 | 0.1988 | 0.6778 | 0.0003 | 0.8756 | 0.4095 |
| 0.7955 | 0.8219 | 0.3813 | 0.1149 | 0.6489 | 0.0002 | 0.6534 | 0.9313 |
| 0.8880 | 0.8885 | 0.2891 | 0.7200 | 0.5242 | 0.0001 | 0.9437 | 0.7806 |
| 0.8566 | 0.6884 | 0.7427 | 0.7622 | 0.6241 | 0.0000 | 0.5855 | 0.4431 |
| 0.9095 | 0.9072 | 0.7596 | 0.1121 | 0.4218 | 0.0000 | 0.6953 | 0.7361 |
| 0.8457 | 0.6463 | 0.6731 | 0.2701 | 0.5199 | 0.0003 | 0.9494 | 0.1078 |
| 0.8453 | 0.7249 | 0.5598 | 0.3468 | 0.6246 | 0.0000 | 0.9781 | 0.0551 |
| 0.8901 | 0.6169 | 0.7800 | 0.4173 | 0.6084 | 0.0001 | 0.8081 | 0.6171 |
| 0.8883 | 0.6349 | 0.9140 | 0.4161 | 0.4861 | 0.0002 | 0.9908 | 0.5416 |
| 0.5891 | 0.8686 | 0.1878 | 0.3771 | 0.5730 | 0.0000 | 0.5980 | 0.9694 |
| 0.9767 | 0.7419 | 0.6312 | 0.5234 | 0.6208 | 0.0001 | 0.9807 | 0.7367 |
| 0.8819 | 0.8242 | 0.7620 | 0.5085 | 0.5998 | 0.0009 | 0.9017 | 0.1560 |
| 0.9689 | 0.7608 | 0.4916 | 0.8493 | 0.6597 | 0.0003 | 0.7385 | 0.6987 |
| 0.9374 | 0.7451 | 0.7940 | 0.6464 | 0.5986 | 0.0000 | 0.8528 | 0.6862 |
| 0.8512 | 0.8238 | 0.6989 | 0.3478 | 0.7022 | 0.0000 | 0.8209 | 0.6218 |
| 0.9488 | 0.7285 | 0.5193 | 0.3102 | 0.4233 | 0.0002 | 0.8902 | 0.0467 |
| 0.7607 | 0.8930 | 0.5969 | 0.2487 | 0.6883 | 0.0001 | 0.9874 | 0.5072 |
| 0.8612 | 0.5033 | 0.8851 | 0.2985 | 0.6186 | 0.0000 | 0.9051 | 0.8365 |
| 0.8661 | 0.9068 | 0.7004 | 0.4906 | 0.5126 | 0.0000 | 0.6797 | 0.2837 |
| 0.9584 | 0.8042 | 0.9085 | 0.1706 | 0.7058 | 0.0000 | 0.8020 | 0.1470 |
| 0.8282 | 0.7414 | 0.3772 | 0.5962 | 0.6157 | 0.0000 | 0.7915 | 0.4168 |
| 0.8080 | 0.4806 | 0.3755 | 0.5792 | 0.7095 | 0.0000 | 0.9208 | 0.9153 |
| 0.8827 | 0.8390 | 0.7363 | 0.4322 | 0.5868 | 0.0002 | 0.9547 | 0.5613 |
| 0.9003 | 0.6074 | 0.9254 | 0.7219 | 0.5336 | 0.0000 | 0.9187 | 0.7110 |
| 0.7717 | 0.4216 | 0.7186 | 0.7085 | 0.5618 | 0.0001 | 0.9906 | 0.7768 |
| 0.8653 | 0.9469 | 0.8303 | 0.6427 | 0.5769 | 0.0000 | 0.9723 | 0.8612 |
| 0.7309 | 0.8329 | 0.6769 | 0.4444 | 0.5761 | 0.0001 | 0.9815 | 0.1824 |
| 0.8133 | 0.8988 | 0.4800 | 0.2895 | 0.5899 | 0.0000 | 0.9815 | 0.7976 |
| 0.8936 | 0.7152 | 0.5950 | 0.7803 | 0.6126 | 0.0000 | 0.9814 | 0.7396 |
| 0.9148 | 0.4461 | 0.5953 | 0.2127 | 0.6978 | 0.0001 | 0.8441 | 0.3009 |
| 0.9275 | 0.9482 | 0.5661 | 0.2924 | 0.5828 | 0.0000 | 0.7495 | 0.6490 |
| 0.6791 | 0.8607 | 0.3370 | 0.5437 | 0.5818 | 0.0000 | 0.8619 | 0.2562 |
| 0.7788 | 0.8872 | 0.5104 | 0.7009 | 0.5262 | 0.0001 | 0.9764 | 0.2308 |
| 0.9081 | 0.8565 | 0.6190 | 0.9002 | 0.6600 | 0.0000 | 0.9739 | 0.6210 |
| 0.9031 | 0.8101 | 0.6756 | 0.7886 | 0.6351 | 0.0000 | 0.5604 | 0.3443 |
| 0.9088 | 0.9131 | 0.9030 | 0.6526 | 0.7004 | 0.0001 | 0.8978 | 0.7916 |
| 0.9121 | 0.7509 | 0.8345 | 0.6139 | 0.6640 | 0.0003 | 0.9085 | 0.4244 |
| 0.9811 | 0.7986 | 0.4625 | 0.0877 | 0.7252 | 0.0000 | 0.7438 | 0.9819 |
| 0.8524 | 0.8802 | 0.7692 | 0.5650 | 0.6305 | 0.0000 | 0.9904 | 0.1062 |
| 0.8026 | 0.9458 | 0.5296 | 0.5765 | 0.6311 | 0.0000 | 0.5763 | 0.1039 |
| 0.7936 | 0.6969 | 0.7149 | 0.5417 | 0.5862 | 0.0000 | 0.9404 | 0.2871 |
| 0.8969 | 0.8723 | 0.7254 | 0.4435 | 0.5442 | 0.0000 | 0.8534 | 0.2496 |
| 0.9145 | 0.7404 | 0.5182 | 0.1573 | 0.5665 | 0.0000 | 0.6989 | 0.7670 |
| 0.7304 | 0.4240 | 0.6651 | 0.7343 | 0.6460 | 0.0005 | 0.8883 | 0.3680 |
| 0.8499 | 0.6778 | 0.5321 | 0.5860 | 0.6003 | 0.0000 | 0.9115 | 1.0000 |
| 0.8107 | 0.5746 | 0.4475 | 0.8809 | 0.5846 | 0.0000 | 0.9522 | 0.7477 |
| 0.9530 | 0.7455 | 0.6706 | 0.2152 | 0.5506 | 0.0000 | 0.8941 | 0.8640 |
| 0.9811 | 0.7049 | 0.7569 | 0.4296 | 0.7359 | 0.0000 | 0.8164 | 0.8569 |
| 0.8282 | 0.7039 | 0.6899 | 0.2621 | 0.6497 | 0.0001 | 0.9799 | 0.1912 |
| 0.8685 | 0.6832 | 0.4910 | 0.1205 | 0.6786 | 0.0000 | 0.8611 | 0.7957 |
| 0.9452 | 0.6914 | 0.1355 | 0.8832 | 0.6140 | 0.0000 | 0.7382 | 0.1247 |
| 0.8891 | 0.7970 | 0.8864 | 0.7121 | 0.6440 | 0.0001 | 0.8486 | 0.1435 |
| 0.8628 | 0.8557 | 0.4464 | 0.4487 | 0.3550 | 0.0000 | 0.6617 | 0.4975 |
| 0.6847 | 0.8768 | 0.5704 | 0.7779 | 0.4518 | 0.0001 | 0.7889 | 0.0538 |
| 0.8160 | 0.5285 | 0.7451 | 0.9605 | 0.6499 | 0.0001 | 0.7973 | 0.9004 |
| 0.7015 | 0.7866 | 0.6060 | 0.1176 | 0.6012 | 0.0000 | 0.9775 | 0.0913 |
| 0.8962 | 0.9211 | 0.7915 | 0.2750 | 0.7283 | 0.0000 | 0.8701 | 0.6017 |
| 0.7408 | 0.8412 | 0.7540 | 0.7592 | 0.5258 | 0.0001 | 0.6898 | 0.6433 |
| 0.7998 | 0.8443 | 0.5343 | 0.2858 | 0.5327 | 0.0000 | 0.9968 | 0.7246 |
| 0.8356 | 0.8412 | 0.5683 | 0.3119 | 0.5914 | 0.0000 | 0.9559 | 0.9975 |
| 0.8477 | 0.7801 | 0.5209 | 0.9154 | 0.6963 | 0.0000 | 0.7167 | 0.9491 |
| 0.7344 | 0.8381 | 0.7348 | 0.2063 | 0.5581 | 0.0000 | 0.9717 | 0.9782 |
| 0.8686 | 0.8368 | 0.7307 | 0.3324 | 0.5977 | 0.0001 | 0.9519 | 0.0875 |
| 0.8200 | 0.5152 | 0.7277 | 0.0821 | 0.5894 | 0.0000 | 0.9561 | 0.0325 |
| 0.8938 | 0.8919 | 0.7762 | 0.3386 | 0.4932 | 0.0001 | 0.9800 | 0.1551 |
| 0.8734 | 0.8803 | 0.9293 | 0.5963 | 0.7255 | 0.0000 | 0.8580 | 0.6905 |
| 0.6591 | 0.9007 | 0.4542 | 0.7206 | 0.5821 | 0.0001 | 0.8649 | 0.9226 |
| 0.9368 | 0.8350 | 0.4827 | 0.3746 | 0.6110 | 0.0001 | 0.8501 | 0.0976 |
| 0.9041 | 0.8029 | 0.5129 | 0.4846 | 0.6182 | 0.0000 | 0.8464 | 0.8915 |
| 0.8124 | 0.8787 | 0.6889 | 0.6279 | 0.5835 | 0.0001 | 0.6100 | 0.9451 |
| 0.8529 | 0.7622 | 0.6015 | 0.3832 | 0.5891 | 0.0000 | 0.9892 | 0.8282 |
| 0.9138 | 0.6880 | 0.5724 | 0.6426 | 0.7034 | 0.0000 | 0.9000 | 0.9257 |
| 0.6823 | 0.8729 | 0.4270 | 0.6149 | 0.7211 | 0.0000 | 0.8255 | 0.1217 |
| 0.7773 | 0.6954 | 0.6682 | 0.8490 | 0.6204 | 0.0000 | 0.8360 | 0.9153 |
| 0.9796 | 0.5894 | 0.6735 | 0.7701 | 0.6290 | 0.0000 | 0.6712 | 0.6244 |
| 0.7655 | 0.8020 | 0.2996 | 0.6398 | 0.7290 | 0.0000 | 0.6556 | 0.7536 |
| 0.8493 | 0.8712 | 0.7570 | 0.3153 | 0.6859 | 0.0000 | 0.8827 | 0.7618 |
| 0.7349 | 0.9371 | 0.2763 | 0.3129 | 0.5076 | 0.0000 | 0.9568 | 0.5302 |
| 0.8280 | 0.7601 | 0.7271 | 0.0425 | 0.4937 | 0.0000 | 0.9843 | 0.6939 |
| 0.9141 | 0.9019 | 0.8922 | 0.0877 | 0.6878 | 0.0002 | 0.9094 | 0.5370 |
| 0.8558 | 0.7342 | 0.6582 | 0.7090 | 0.5759 | 0.0000 | 0.9905 | 0.2880 |
| 0.8323 | 0.8098 | 0.6159 | 0.7092 | 0.4719 | 0.0001 | 0.8332 | 0.4014 |
| 0.9112 | 0.4960 | 0.7879 | 0.8688 | 0.6448 | 0.0000 | 0.5576 | 0.3525 |
| 0.9400 | 0.7532 | 0.8286 | 0.4820 | 0.6995 | 0.0000 | 0.9087 | 0.9761 |
| 0.8445 | 0.6109 | 0.5063 | 0.3412 | 0.5388 | 0.0000 | 0.9992 | 0.8313 |
| 0.7862 | 0.8018 | 0.7787 | 0.6946 | 0.5632 | 0.0000 | 0.9381 | 0.7039 |
| 0.9087 | 0.8532 | 0.7890 | 0.9865 | 0.6511 | 0.0001 | 0.8940 | 0.1467 |
| 0.9222 | 0.5808 | 0.8795 | 0.2853 | 0.7204 | 0.0000 | 0.8595 | 0.8354 |
| 0.5958 | 0.7108 | 0.6615 | 0.9086 | 0.6875 | 0.0000 | 0.7952 | 0.8135 |
| 0.8688 | 0.7178 | 0.4262 | 0.8309 | 0.6432 | 0.0000 | 0.8940 | 0.0986 |
| 0.9378 | 0.6855 | 0.4334 | 0.7848 | 0.6220 | 0.0001 | 0.8917 | 0.5322 |
| 0.8393 | 0.7849 | 0.6135 | 0.6264 | 0.5498 | 0.0002 | 0.9147 | 0.2098 |
| 0.7993 | 0.7493 | 0.6278 | 0.7945 | 0.6068 | 0.0000 | 0.9089 | 0.7555 |
| 0.6895 | 0.8827 | 0.6245 | 0.5696 | 0.4567 | 0.0000 | 0.8812 | 0.1238 |
| 0.9311 | 0.9071 | 0.9375 | 0.3701 | 0.5791 | 0.0000 | 0.9931 | 0.2271 |
| 0.8853 | 0.5901 | 0.7807 | 0.6472 | 0.5224 | 0.0000 | 0.8016 | 0.6924 |
| 0.9211 | 0.5338 | 0.9045 | 0.1731 | 0.6599 | 0.0005 | 0.8756 | 0.2381 |
| 0.9189 | 0.9058 | 0.6033 | 0.4776 | 0.4806 | 0.0000 | 0.9398 | 0.5310 |
| 0.9571 | 0.7843 | 0.7608 | 0.0569 | 0.7143 | 0.0000 | 0.8079 | 0.8993 |
| 0.7320 | 0.8187 | 0.4132 | 0.1612 | 0.4864 | 0.0003 | 0.8594 | 0.8005 |
| 0.8628 | 0.8282 | 0.3843 | 0.1334 | 0.7237 | 0.0005 | 0.7711 | 0.4556 |
| 0.8424 | 0.7820 | 0.4851 | 0.5230 | 0.5902 | 0.0000 | 0.7813 | 0.0939 |
| 0.7528 | 0.7999 | 0.5523 | 0.7859 | 0.5783 | 0.0006 | 0.7939 | 0.6439 |
| 0.9269 | 0.5708 | 0.8656 | 0.6522 | 0.5871 | 0.0000 | 0.6341 | 0.5674 |
| 0.9437 | 0.8355 | 0.8201 | 0.6189 | 0.7418 | 0.0000 | 0.8027 | 0.5086 |
| 0.9145 | 0.7825 | 0.7177 | 0.8639 | 0.6438 | 0.0000 | 0.7847 | 0.8704 |
| 0.7793 | 0.8719 | 0.4677 | 0.8158 | 0.4985 | 0.0000 | 0.6418 | 0.3760 |
| 0.8587 | 0.7913 | 0.8439 | 0.4875 | 0.6241 | 0.0001 | 0.9399 | 0.8492 |
| 0.9211 | 0.9434 | 0.6736 | 0.4147 | 0.4618 | 0.0001 | 0.9181 | 0.1100 |
| 0.7924 | 0.8475 | 0.5815 | 0.3142 | 0.5756 | 0.0002 | 0.9491 | 0.8740 |
| 0.9448 | 0.8760 | 0.7077 | 0.4483 | 0.7551 | 0.0004 | 0.7666 | 0.0517 |
| 0.9493 | 0.7879 | 0.6666 | 0.1367 | 0.6996 | 0.0000 | 0.7622 | 0.8988 |
| 0.9362 | 0.6982 | 0.7643 | 0.6638 | 0.5505 | 0.0002 | 0.9406 | 0.8437 |
| 0.8563 | 0.4065 | 0.4975 | 0.5676 | 0.7669 | 0.0001 | 0.8919 | 0.3961 |
| 0.8821 | 0.8480 | 0.4214 | 0.7630 | 0.5454 | 0.0000 | 0.9339 | 0.9140 |
| 0.8909 | 0.6577 | 0.6279 | 0.3132 | 0.5284 | 0.0001 | 0.9140 | 0.3556 |
| 0.9635 | 0.5664 | 0.7039 | 0.3465 | 0.7229 | 0.0000 | 0.7142 | 0.9032 |
| 0.9214 | 0.8755 | 0.6663 | 0.6914 | 0.8527 | 0.0000 | 0.3443 | 0.9229 |
| 0.9672 | 0.7731 | 0.2508 | 0.4833 | 0.6151 | 0.0003 | 0.8880 | 0.7715 |
| 0.5258 | 0.8689 | 0.5089 | 0.2685 | 0.4758 | 0.0001 | 0.9318 | 0.8786 |
| 0.8472 | 0.8149 | 0.4629 | 0.2935 | 0.6183 | 0.0001 | 0.9188 | 0.8775 |
| 0.9175 | 0.9232 | 0.6697 | 0.4866 | 0.4545 | 0.0000 | 0.7849 | 0.5128 |
| 0.9216 | 0.6884 | 0.9324 | 0.3853 | 0.7130 | 0.0000 | 0.8530 | 0.5370 |
| 0.9821 | 0.8377 | 0.8336 | 0.3274 | 0.6218 | 0.0000 | 0.9140 | 0.7278 |
| 0.9268 | 0.6704 | 0.7000 | 0.1591 | 0.5256 | 0.0001 | 0.8709 | 0.9692 |
| 0.9485 | 0.8869 | 0.7191 | 0.0498 | 0.4324 | 0.0000 | 0.8763 | 0.7742 |
| 0.5443 | 0.6981 | 0.4967 | 0.8802 | 0.6047 | 0.0000 | 0.9319 | 0.9630 |
| 0.5299 | 0.8746 | 0.5000 | 0.5695 | 0.7106 | 0.0001 | 0.4684 | 0.0344 |
| 0.7717 | 0.8328 | 0.4451 | 0.6480 | 0.5376 | 0.0000 | 0.8576 | 0.0037 |
| 0.9457 | 0.8873 | 0.8164 | 0.9801 | 0.7593 | 0.0008 | 0.6038 | 0.8127 |
| 0.9283 | 0.7913 | 0.8500 | 0.9711 | 0.7500 | 0.0000 | 0.6446 | 0.8340 |
| 0.8842 | 0.4609 | 0.5968 | 0.3135 | 0.4802 | 0.0001 | 0.8639 | 0.9903 |
| 0.7559 | 0.6760 | 0.6238 | 0.2374 | 0.6457 | 0.0001 | 0.6886 | 0.3122 |
| 0.8860 | 0.8828 | 0.6218 | 0.9229 | 0.5415 | 0.0003 | 0.7047 | 0.0674 |
| 0.9666 | 0.5586 | 0.6984 | 0.8792 | 0.4863 | 0.0001 | 0.8975 | 0.5899 |
| 0.9675 | 0.8490 | 0.8137 | 0.0866 | 0.6646 | 0.0001 | 0.8138 | 0.1799 |
| 0.8861 | 0.6183 | 0.8538 | 0.6442 | 0.6438 | 0.0001 | 0.9753 | 0.0258 |
| 0.8719 | 0.8630 | 0.5634 | 0.7399 | 0.6751 | 0.0000 | 0.8562 | 0.0992 |
| 0.9595 | 0.5519 | 0.7457 | 0.6456 | 0.5487 | 0.0000 | 0.8333 | 0.6328 |
| 0.9932 | 0.6806 | 0.5158 | 0.2099 | 0.6042 | 0.0001 | 0.8611 | 0.0984 |
| 0.8030 | 0.8894 | 0.7054 | 0.3233 | 0.5865 | 0.0002 | 0.9608 | 0.0736 |
| 0.6854 | 0.9117 | 0.6404 | 0.2875 | 0.4698 | 0.0000 | 0.9725 | 0.2106 |
| 0.9868 | 0.7832 | 0.8280 | 0.5208 | 0.6530 | 0.0000 | 0.6549 | 0.4411 |
| 0.8715 | 0.5164 | 0.7120 | 0.3223 | 0.5960 | 0.0001 | 0.8495 | 0.0087 |
| 0.7771 | 0.8400 | 0.4902 | 0.5355 | 0.6429 | 0.0004 | 0.8412 | 0.5957 |
| 0.8606 | 0.6150 | 0.7503 | 0.1583 | 0.5434 | 0.0001 | 0.6634 | 0.5453 |
| 0.8783 | 0.8142 | 0.7442 | 0.3184 | 0.6410 | 0.0000 | 0.9012 | 0.9830 |
| 0.8624 | 0.9277 | 0.7179 | 0.8103 | 0.5092 | 0.0000 | 0.9903 | 0.9650 |
| 0.9584 | 0.6816 | 0.8197 | 0.3028 | 0.6217 | 0.0000 | 0.8337 | 0.9085 |
| 0.9546 | 0.8183 | 0.9097 | 0.1786 | 0.6169 | 0.0000 | 0.8507 | 0.3122 |
| 0.7940 | 0.7605 | 0.5514 | 0.3142 | 0.5925 | 0.0001 | 0.9984 | 0.8406 |
| 0.6550 | 0.5904 | 0.7017 | 0.8935 | 0.4274 | 0.0000 | 0.9018 | 0.2682 |
| 0.9265 | 0.7982 | 0.8154 | 0.7682 | 0.4721 | 0.0001 | 0.9689 | 0.0276 |
| 0.6824 | 0.7910 | 0.6963 | 0.4296 | 0.6807 | 0.0000 | 0.7864 | 0.9979 |
| 0.8796 | 0.6090 | 0.5000 | 0.6296 | 0.6089 | 0.0000 | 0.9922 | 0.7191 |
| 0.8563 | 0.7983 | 0.8874 | 0.8546 | 0.7268 | 0.0000 | 0.9441 | 0.9362 |
| 0.8988 | 0.4811 | 0.6169 | 0.8097 | 0.6536 | 0.0001 | 0.6416 | 0.8357 |
| 0.9529 | 0.8047 | 0.5319 | 0.4389 | 0.4421 | 0.0000 | 0.9998 | 0.5381 |
| 0.7421 | 0.8275 | 0.2600 | 0.8473 | 0.6721 | 0.0001 | 0.9157 | 0.8246 |
| 0.7395 | 0.8362 | 0.7389 | 0.6610 | 0.3660 | 0.0000 | 0.6871 | 0.9268 |
| 0.9707 | 0.6154 | 0.7218 | 0.4491 | 0.6774 | 0.0003 | 0.9557 | 0.0615 |
| 0.9582 | 0.7463 | 0.8911 | 0.6920 | 0.6368 | 0.0004 | 0.6624 | 0.2963 |
| 0.8751 | 0.7825 | 0.8157 | 0.5978 | 0.5991 | 0.0001 | 0.9269 | 0.7503 |
| 0.9050 | 0.7909 | 0.6362 | 0.8512 | 0.7158 | 0.0001 | 0.8674 | 0.6917 |
| 0.8813 | 0.8737 | 0.2587 | 0.4873 | 0.7013 | 0.0000 | 0.8917 | 0.0718 |
| 0.7802 | 0.7807 | 0.5714 | 0.5087 | 0.6777 | 0.0001 | 0.5428 | 0.0445 |
| 0.9575 | 0.9042 | 0.7286 | 0.0677 | 0.6995 | 0.0000 | 0.7927 | 0.1669 |
| 0.9636 | 0.6325 | 0.8507 | 0.3775 | 0.6177 | 0.0000 | 0.9514 | 0.7745 |
| 0.4927 | 0.6686 | 0.6517 | 0.1575 | 0.5890 | 0.0001 | 0.9927 | 0.9924 |
| 0.9341 | 0.7209 | 0.7195 | 0.6463 | 0.6738 | 0.0000 | 0.8913 | 0.9214 |
| 0.9071 | 0.7330 | 0.4839 | 0.4205 | 0.8446 | 0.0000 | 0.4691 | 0.9577 |
| 0.8930 | 0.8804 | 0.4262 | 0.4809 | 0.5866 | 0.0001 | 0.9941 | 0.1176 |
| 0.8228 | 0.7032 | 0.8236 | 0.6142 | 0.5868 | 0.0001 | 0.9652 | 0.5512 |
| 0.7350 | 0.8202 | 0.8231 | 0.7414 | 0.5969 | 0.0007 | 0.9402 | 0.2761 |
| 0.7208 | 0.5817 | 0.8436 | 0.9441 | 0.5253 | 0.0000 | 0.8649 | 0.9252 |
| 0.8918 | 0.7716 | 0.8144 | 0.0504 | 0.5576 | 0.0000 | 0.9972 | 0.8009 |
| 0.8046 | 0.8155 | 0.7076 | 0.0229 | 0.6293 | 0.0000 | 0.8328 | 0.8469 |
| 0.8965 | 0.5716 | 0.8519 | 0.8638 | 0.6446 | 0.0001 | 0.4906 | 0.8231 |
| 0.5775 | 0.7676 | 0.6667 | 0.2012 | 0.5418 | 0.0001 | 0.9635 | 0.5747 |
| 0.9054 | 0.9568 | 0.7530 | 0.4442 | 0.6778 | 0.0001 | 0.9672 | 0.8523 |
| 0.8516 | 0.6741 | 0.5624 | 0.9219 | 0.7224 | 0.0000 | 0.7850 | 0.3196 |
| 0.9184 | 0.5130 | 0.6712 | 0.9485 | 0.5872 | 0.0000 | 0.9835 | 0.9998 |
| 0.9632 | 0.8548 | 0.7763 | 0.1924 | 0.7233 | 0.0000 | 0.8939 | 0.2976 |
| 0.9259 | 0.9925 | 0.7632 | 0.2719 | 0.7080 | 0.0001 | 0.7114 | 0.7216 |
| 0.8122 | 0.8400 | 0.6794 | 0.2307 | 0.5825 | 0.0000 | 0.9891 | 0.9021 |
| 0.8944 | 0.6586 | 0.8692 | 0.5460 | 0.5845 | 0.0000 | 0.8598 | 0.1533 |
| 0.9101 | 0.8200 | 0.9436 | 0.3334 | 0.6593 | 0.0000 | 0.5019 | 0.1892 |
| 0.9343 | 0.7303 | 0.8252 | 0.4582 | 0.5375 | 0.0000 | 0.9917 | 0.3141 |
| 0.7756 | 0.7334 | 0.6753 | 0.4121 | 0.5570 | 0.0000 | 0.7474 | 0.3721 |
| 0.8910 | 0.9057 | 0.7637 | 0.1398 | 0.7162 | 0.0000 | 0.6734 | 0.6823 |

The histograms of 1000 replicated Z-residual test p-values for the wbc
and lwbc models. The red vertical lines in these histograms show the
upper bound summaries of these replicated p-values, p\_{min}. These
histograms show that the Z-SW, Z-SF, and Z-AOV with LP tests for both
models give a large proportion of p-values greater than 0.05, and the
large p-values result in large p\_{min} values. In contrast, the
replicated Z-AOV with log(wbc) p-values for the lwbc model are almost
all smaller than 0.001. The consistently small Z-AOV with log(wbc)
p-values further confirm that the log transformation of wbc is
inappropriate for modelling the survival time.

[TABLE]

Code

``` r
par(mfrow = c(4, 2), mar = c(4, 4, 2, 2))

hist(
  sw.wbc,
  main  = "Replicated Z-SW P-values for wbc Model",
  breaks = 20,
  xlab  = "Z-SW P-values for wbc Model"
)
abline(v = pmin.sw.LeukSurv.wbc, col = "red")

hist(
  sw.lwbc,
  main  = "Replicated Z-SW P-values for lwbc Model",
  breaks = 20,
  xlab  = "Z-SW P-values for lwbc Model"
)
abline(v = pmin.sw.LeukSurv.lwbc, col = "red")

hist(
  sf.wbc,
  main  = "Replicated Z-SF P-values for wbc Model",
  breaks = 20,
  xlab  = "Z-SF P-values for wbc Model"
)
abline(v = pmin.sf.LeukSurv.wbc, col = "red")

hist(
  sf.lwbc,
  main  = "Replicated Z-SF P-values for lwbc Model",
  breaks = 20,
  xlab  = "Z-SF P-values for lwbc Model"
)
abline(v = pmin.sf.LeukSurv.lwbc, col = "red")

hist(
  aov.wbc.lp,
  main  = "Replicated Z-AOV with LP P-values for wbc Model",
  breaks = 20,
  xlab  = "Z-AOV with LP P-values for wbc Model"
)
abline(v = pmin.aov.lp.LeukSurv.wbc, col = "red")

hist(
  aov.lwbc.lp,
  main  = "Replicated Z-AOV with LP P-values for lwbc Model",
  breaks = 20,
  xlab  = "Z-AOV with LP P-values for lwbc Model"
)
abline(v = pmin.aov.lp.LeukSurv.lwbc, col = "red")

hist(
  aov.wbc,
  main  = "Replicated Z-AOV with wbc P-values for wbc Model",
  breaks = 20,
  xlab  = "Z-AOV with wbc P-values for wbc Model"
)
abline(v = pmin.aov.wbc.LeukSurv, col = "red")

hist(
  aov.lwbc,
  main  = "Replicated Z-AOV with wbc P-values for wbc Model",
  breaks = 20,
  xlab  = "Z-AOV with wbc P-values for lwbc Model"
)
abline(v = pmin.aov.lwbc.LeukSurv, col = "red")
```

![](zresidual_demo_files/figure-html/fig-zresid-hist-pvalues-1.png)

Figure 3: Figure 5: The histograms of 1000 replicated Z-SW, Z-SF,
Z-AOV-LP and Z-AOV-log(wbc) p-values for the wbc model (left panels) and
the lwbc model (right panels) fitted with the survival times of acute
myeloid leukemia patients. The vertical red lines indicate p\_{min} for
1000 replicated p-values. Note that the upper limit of the x-axis for
Z-AOV-log(wbc) p-values for the lwbc model is 0.005, not 1 for others.

## 5 Other residual calculation

### 5.1 censored Z-residuals

The normality of censored Z-residuals is tested by an extended SF method
for censored observations, which is implemented with gofTestCensored in
the R package EnvStats.

Code

``` r
censored.Zresid.LeukSurv.wbc<-surv_residuals(fit.object = fit_LeukSurv_wbc,
                                      data=LeukSurv,
                                      residual.type="censored Z-residual")

censored.Zresid.LeukSurv.logwbc<-surv_residuals(fit.object = fit_LeukSurv_logwbc,data= LeukSurv,residual.type="censored Z-residual")

gof.censored.zresidual(censored.Zresidual=censored.Zresid.LeukSurv.wbc)
```

    [1] 0.5702324

Code

``` r
gof.censored.zresidual(censored.Zresidual=censored.Zresid.LeukSurv.logwbc)
```

    [1] 0.07535993

### 5.2 Cox-Snell residual

The overall GOF tests and graphical checking with Cox-Snell residuals
show that both the wbc and lwbc models provide adequate fits to the
dataset. The estimated CHFs of the CS residuals of both of the wbc and
lwbc models align closely along the 45^{\circ} diagonal line.

Code

``` r
##unmodified CS residuals
ucs.LeukSurv.wbc<-surv_residuals(fit.object = fit_LeukSurv_wbc,data= LeukSurv,residual.type = "Cox-Snell" )
ucs.LeukSurv.logwbc<-surv_residuals(fit.object = fit_LeukSurv_logwbc,data= LeukSurv,residual.type = "Cox-Snell" )
```

Code

``` r
##unmodified CS residuals
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
plot.cs.residual(ucs.LeukSurv.wbc,main.title = "CS Residuals of wbc model")
plot.cs.residual(ucs.LeukSurv.logwbc,main.title = "CS Residuals of lwbc model")
```

![](zresidual_demo_files/figure-html/unnamed-chunk-6-1.png)

### 5.3 Martingale residual

The martingale residuals are mostly within the interval (-3, 1) for
those two models. In the scatterplots of martingale residuals under the
wbc model, the LOWESS curves have a slight upward slope on the left,
while under the lwbc model, they display a pronounced downward curve.
Both of these lines demonstrate noticeable non-horizontal trends.

Code

``` r
martg.LeukSurv.wbc<-surv_residuals(fit.object = fit_LeukSurv_wbc,data= LeukSurv,residual.type = "martingale")
martg.LeukSurv.logwbc<-surv_residuals(fit.object = fit_LeukSurv_logwbc,data= LeukSurv,residual.type = "martingale" )
```

Code

``` r
par(mfrow = c(1,2))
plot.martg.resid(martg.LeukSurv.wbc,x_axis_var="wbc",main.title = "Martingale Residuals of wbc Model")
plot.martg.resid(martg.LeukSurv.logwbc,x_axis_var="logwbc",main.title = "Martingale Residuals of lwbc Model")
```

![](zresidual_demo_files/figure-html/unnamed-chunk-8-1.png)

### 5.4 Deviance residual

The deviance residuals are more symmetrically distributed than
martingale residuals and they are mostly within the interval (-3, 3). In
both models, the scatterplots of deviance residuals exhibit strikingly
non-horizontal trends in their LOWESS curves.

Code

``` r
#Deviance residuals
dev.LeukSurv.wbc<-surv_residuals(fit.object = fit_LeukSurv_wbc,data= LeukSurv,residual.type = "deviance" )
dev.LeukSurv.logwbc<-surv_residuals(fit.object = fit_LeukSurv_logwbc,data= LeukSurv,residual.type = "deviance" )
```

Code

``` r
par(mfrow = c(1,2))
plot.dev.resid(dev.LeukSurv.wbc,x_axis_var="wbc",main.title = "Deviance Residuals of wbc Model")
plot.dev.resid(dev.LeukSurv.logwbc,x_axis_var="logwbc",main.title = "Deviance Residuals of lwbc Model")
```

![](zresidual_demo_files/figure-html/unnamed-chunk-10-1.png)

## 6 References

Wu, T., Li, L., & Feng, C. (2024). Z-residual diagnostic tool for
assessing covariate functional form in shared frailty models. Journal of
Applied Statistics, 52(1), 28–58.
https://doi.org/10.1080/02664763.2024.2355551
