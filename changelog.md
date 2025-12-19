## Change Logs

### December 19, 2025

* In the vignettes, you can generate a large number like 1000 replicated Z-residuals so that the histograms look better, and the upper-bound p-values have more sample size. using the saved rds to save the pre-computed objects to inst/extdata to save the re-production time. See this: https://gemini.google.com/share/5362ba39f503 

* Use the upper-bound p-values for the median to summarize the replicated Z-residuals for all of them, even for the past works on survival Z-residuals. This is a more theoretically grounded method than p-min. Don’t use the min p-value, which confuses people. Use the name in the recent paper on this.

* Sanitize the function names so that the names follow the same patterns, for example using log_pred.bern.brms, phurdle, ptruncnb, etc. The functions like pdf and cdf needs to be renamed to be consistent with conventional names dnorm.pnorm.  If the function is applied to a class name, register an S3 method. You can use either of these two methods:
  - Rename (refactor) all the messy function names. But need to change all of them in the packages. In RStudio, you can use edit-> find in files to find all the function names in the whole directory. This may be the easiest way to change all these function names. I would recommend this permanent change for future developers.
  - Using S3method export to create an alias for users. This may be easier without modifying the internal code but may cause confusion to future developers.

*
  - Write a method to compute log_pred given fitting class: log_pred.brms, log_pred.coxph,log_pred.survreg, 

  - Write a method Zresidual for class log_pred (log_cdf, log_pmf, log_like), all **vectors**

* Improve vignettes like: $$ $${#eq-rsp}, more verbal descriptions of the experiments and results.

