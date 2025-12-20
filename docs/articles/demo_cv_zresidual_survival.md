# Demo of Cross-validatory Z-Residual for Diagnosing Shared Frailty Models

## 1 Installing Zresidual and Other packages

### 1.1 Installing Z-residua from the source

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
cross-validatory (CV) Z-residuals based on the output of the coxph
function from the survival package in R. It also serves as a
demonstration of how to use cross-validatory Z-residuals to identify
outliers in semi-parametric shared frailty models. To fully understand
the detailed definitions and the example data analysis results, please
refer to the original paper titled “Cross-validatory Z-Residual for
Diagnosing Shared Frailty Models”.

## 3 Definition of Cross-validatory Z-residual

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
given x\_{i} under the true model. Therefore, the RSPs can be
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

In our study, we employ both leave-one-out cross-validation (LOOCV) and
10-fold cross-validation (10-fold CV) techniques to compute
cross-validatory Z-residuals. In LOOCV Z-residual, one observation,
t\_{ij}^{test}, is excluded from the dataset with n observations. The
remaining observations, acting as the training dataset, are used for
parameter estimation in the shared frailty model. Fitting the model to
the training dataset produces the estimated regression coefficients,
\hat{\beta'}, and frailty effects, \hat{z_i}. The Breslow estimator
helps estimate the cumulative baseline hazard (\hat{H_0}). The survival
function \hat{S}{ij} (y{ij}) for the test observation y\_{ij}^{test} is
computed using: \begin{equation} \hat{S}\_{ij}(y\_{ij}^{test}) = \exp
\\- \hat{z_i} \exp(\hat{\beta'} x\_{ij}) \hat{H}\_0(y\_{ij}^{test}) \\.
\end{equation} Subsequently, the RSP for the observed t\_{ij} is defined
as: \begin{equation} \hat{S}\_{ij}^{R}(t\_{ij}^{test}, d\_{ij},
U\_{ij})= \left\\ \begin{array}{rl} \hat{S}\_{ij}(t\_{ij}^{test}), &
\text{if \$t\_{ij}^{test}\$ is uncensored, i.e., \$d\_{ij}=1\$,}\\
U\_{ij}\\\hat{S}\_{ij}(t\_{ij}^{test}), & \text{if \$t\_{ij}^{test}\$ is
censored, i.e., \$d\_{ij}=0\$.} \end{array} \right. \end{equation}
Resulting in the Z-residual for t\_{ij}^{test}: \begin{equation}
\hat{z}\_{ij}(t\_{ij}^{test}, d\_{ij}, U\_{ij})=-\Phi^{-1}
(\hat{S}\_{ij}^R(t\_{ij}^{test}, d\_{ij}, U\_{ij})). \end{equation} By
repeating these steps for each observation (n times), a LOOCV predictive
Z-residual is calculated for each observation.

For cluster-based or categorical covariate values, specific
considerations are employed during the LOOCV and k-fold CV methods.
Clusters with only one observation cannot be included in the training
dataset, and similar requirements are imposed on categorical covariates.
As such, cross-validatory Z-residuals for these observations are
designated as NA in the implementation.

## 4 Examples for Illustration and Demonstration

### 4.1 Load the real Dataset

This example demonstrates the practical application of cross-validatory
Z-residuals in identifying outliers within a study on kidney infections.
The dataset comprises records from 38 kidney patients using a portable
dialysis machine. It documents the times of the first and second
recurrences of kidney infections for these patients. Each patient’s
survival time is defined as the duration until infection from catheter
insertion. The patient records are considered as clusters due to shared
frailty, signifying the common effect across patients. Instances where
the catheter is removed for reasons other than infection are treated as
censored observations, accounting for 24% of the dataset. The dataset
encompasses 38 patient clusters, with each patient having exactly two
observations, resulting in a total sample size of 76. This dataset is
frequently employed to exemplify shared frailty models.

Code

``` r
data_path <- system.file("extdata", "kidney.rda", package = "Zresidual")
load(data_path)
kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
kidney$sex<-as.factor(kidney$sex)
kidney$id<-as.factor(kidney$id)
```

### 4.2 Fitting Models

We fit a shared gamma frailty model with three covariates: covariates:
age in years, gender (male or female), and four disease types (0=GN,
1=AN, 2=PKD, 3=Other).

Code

``` r
fit_kidney <-coxph(Surv(time, status) ~ age + sex + disease+frailty(id, distribution="gamma"), data= kidney)
```

### 4.3 Z-residual and LOOCV Z-residual calculation

We computed Z-residuals using the No-CV and LOOCV methods for the kidney
infection dataset. Given the similarity in performance between the
10-fold CV and LOOCV Z-residual methods demonstrated in the simulation
studies and the manageable computational load, we focused on the LOOCV
method.

Code

``` r
Zresid.kidney<-Zresidual(fit_kidney,nrep=10)
CVZresid.kidney<-CV.Zresidual(fit_kidney,nrep=10,nfolds = nrow(kidney))
CVZresid.kidney_new <- CVZresid.kidney [-57, ]
```

### 4.4 Inspecting the Normality of Z-Residuals for Checking Overall GOF

The first and second columns of Figure 1 display scatterplots against
the index and QQ plots of the Z-residuals calculated with the No-CV and
LOOCV methods. The No-CV Z-residuals predominantly fall within the range
of -3 to 3, displaying alignment with the 45^\circ straight line in the
QQ plot. The QQ plot of the No-CV Z-residuals indicates a SW p-value of
around 0.70, signifying a well-fitted model to the dataset. Thus, the
diagnostic results using No-CV Z-residuals suggest the suitability of
the shared frailty model for the dataset without identifying any
outliers. However, analysis of the scatterplot of LOOCV Z-residuals
reveals that the Z-residuals of cases labeled 20 and 42 exceed 3. These
instances are considered outliers for the shared frailty model. The QQ
plot of LOOCV Z-residuals displays a noticeable deviation from the
45^\circ straight line, attributed to the considerable Z-residuals of
the two identified outliers. The SW p-value of LOOCV Z-residuals is
notably small, measuring less than 0.01, as evident in the QQ plot. In
summary, the diagnosis results with LOOCV Z-residuals suggest that the
fitted shared frailty model is inadequate for this dataset, and two
cases exhibit excessive Z-residuals, categorized as outliers for this
model.

Code

``` r
for (i in 1:10) {
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plot.zresid(Zresid.kidney,x_axis_var="index", main.title = "Z-residual Scatterplot",
     outlier.return = TRUE,irep=i)
plot.zresid(CVZresid.kidney_new,x_axis_var="index", main.title = "LOOCV Z-residual Scatterplot",
     outlier.return = TRUE,irep=i)
qqnorm.zresid(Zresid.kidney, main.title = "Z-residual QQ plot",irep=i)
qqnorm.zresid(CVZresid.kidney_new, main.title = "LOOCV Z-residual QQ plot",irep=i)
}
```

![](demo_cv_zresidual_survival_files/figure-html/qqplot-zresid-.gif)

Figure 1: Scatterplots and QQ plots of No-CV and LOOCV Z-residuals of
the fitted shared frailty models based on the original kidney infection
dataset.

### 4.5 Diagnostic Tests with Z-residuals

The Shapiro-Wilk (SW) or Shapiro-Francia (SF) normality tests applied to
Z-residuals can be used to numerically test the overall GOF of the
model.

Code

``` r
library(gt)
library(tibble)
library(dplyr)
sw.kidney<-sw.test.zresid(Zresid.kidney)
sw.kidney.cv<-sw.test.zresid(CVZresid.kidney)
sf.kidney<-sf.test.zresid(Zresid.kidney)
sf.kidney.cv<-sf.test.zresid(CVZresid.kidney)
gof_tests<-data.frame(sw.kidney,sw.kidney.cv,sf.kidney,sf.kidney.cv)
gof_tests_table <- gof_tests %>%
  rownames_to_column(var = "Metric") %>%
  gt() %>%
  tab_spanner(
    label = "Shapiro-Wilk",
    columns = c(sw.kidney, sw.kidney.cv)
  ) %>%
  tab_spanner(
    label = "Shapiro-Francia",
    columns = c(sf.kidney, sf.kidney.cv)
  ) %>%
  fmt_number(
    columns = -Metric,
    decimals = 4
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_options(
    column_labels.font.weight = "bold"
  )

gof_tests_table
```

[TABLE]

There exists randomness in the Z-residuals of censored observations,
meaning that different sets of Z-residuals can be generated for the same
dataset using distinct random numbers. Thus, to test the robustness of
the previously conducted diagnosis, we replicated a large number of
realizations of Z-residuals. The Figure 2 exhibits histograms of 1000 SW
test p-values, each derived from a set of No-CV or LOOCV Z-residuals.
More than 95% of the SW p-values for No-CV Z-residuals surpass 0.05,
whereas 100% of the SW p-values for LOOCV Z-residuals fall below 0.05.
This consistency across numerous replications confirms that the
evaluation of the misspecification of the shared frailty model is not
incidental to a specific set of LOOCV Z-residuals but a recurring
conclusion supported by extensive Z-residual replications.

[TABLE]

Code

``` r
par(mfrow = c(2,2),mar=c(4,4,2,2))
hist(sw.kidney,main="Replicated Z-SW P-values, No CV",breaks=20,
     xlab="Z-SW P-values")
abline(v=pmin.sw.kidney,col="red")
hist(sw.kidney.cv,main="Replicated Z-SW P-values, LOOCV",breaks=20,
     xlab="Z-SW P-values")
abline(v=pmin.sw.kidney.cv,col="red")

hist(sf.kidney,main="Replicated Z-SF P-values, No CV",breaks=20,
     xlab="Z-SF P-values")
abline(v=pmin.sf.kidney,col="red")
hist(sf.kidney.cv,main="Replicated Z-SF P-values, LOOCV",breaks=20,
     xlab="Z-SF P-values")
abline(v=pmin.sf.kidney.cv,col="red")
```

![](demo_cv_zresidual_survival_files/figure-html/unnamed-chunk-5-1.png)

## 5 References

Wu, T., Feng, C., & Li, L. (2024). Cross-Validatory Z-Residual for
Diagnosing Shared Frailty Models. The American Statistician, 79(2),
198–211. https://doi.org/10.1080/00031305.2024.2421370
