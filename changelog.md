# Change Logs for Zresidual Package


### Planned Tasks on December 19, 2025

- Improve vignettes like:

  - `$$Z_i =-\Phi^{-1} (RSP_i) $${#eq-rsp}`,
  - improve table, figure format and captions
  - more verbal descriptions of the experiments and results.

- In the vignettes, you can generate a large number like 1000 replicated
  Z-residuals so that the histograms look better, and the upper-bound
  p-values have more sample size. using the saved rds to save the
  pre-computed objects to inst/extdata to save the re-production time.
  See this: <https://gemini.google.com/share/5362ba39f503>

- Use the upper-bound p-values for the median to summarize the
  replicated Z-residuals for all of them, even for the past works on
  survival Z-residuals. This is a more theoretically grounded method
  than p-min. Don’t use the min p-value, which confuses people. Use the
  name in the recent paper on this.

- Sanitize the function names so that the names follow the generic.class
  convention, for example using log_pred.bern.brms, phurdle, ptruncnb,
  etc. The functions like pdf and cdf needs to be renamed to be
  consistent with conventional names dnorm.pnorm. If the function is
  applied to a class name, register an S3 method. Rename (refactor) all
  the messy function names. But need to change all of them in the
  packages. In RStudio, you can use edit-\> find in files to find all
  the function names in the whole directory. This may be the easiest way
  to change all these function names. I would recommend this permanent
  change for future developers. This step will need to be combined with
  future work in “refactoring Zresidual” discussed below.

- Refactor Zresidual method and the subsequent worklows to be extensible
  for custom models

  - The structure of the workflow is as follows:

    - Frequentist models:

      model -\> log_pred_family_pkg -\> log_pred_pkg -\> Zresidual -\>
      Diangnose

    - Bayesian models:

      model -\> log_mcmcpred_family_pkg -\> log_mcmcpred_pkg -\>
      log_pred_pkg -\> Zresidual -\> Diagnose

  - **The steps from log_pred_pkg -\> Zresidual -\> Diagnose is generic.
    What’s difference is computing the obj in log_pred_pkg for different
    family. The naming of the functions should be standardized so that
    the upper-level function/method Zresidual will cascadely search the
    functions for computing log_pred for each family and each pkg.**

  - To refactor Zresidual to use the above workflow, we will need to do
    these:

    - Rewrite and rename the lowest level functions to compute log_pred
      given pkg and family, eg,
      log_pred_bern_brms,log_pred_hurdleNB_brms,
      log_pred_coxph,log_pred_weibull_survreg. The family name and pkg
      names and be retrieve from the model object.

    - Rewrite the method Zresidual to take input from class `log_pred`:
      (`log_cdf`, `log_pmf`, both **vectors**, and optional `data`,
      which includes `lp`, and `covariates`, which is needed for plot
      and stat tests), .

    - To generally support naive cross-validation or external residuals,
      in the step from model to log_pred_family_pkg or
      log_mcmcpred_family_pkg, we should allow a testdata arugment,
      which could be an external `testdata`. If `testdata` is null, it
      will be the `data` from model obj, ie, training data

    - It is good to define %\>% operator for the above workflow. We also
      provide short-cut function like Zredidual(model, …) to do these
      steps.

    - If the obj to Zresidual is a specific model (e.g, brms for hurdle
      etc.), find or call lower level functions directly, such as
      log_pred.brms, to compute Z-residuals. Before finishing the whole
      refactoring process, we can leave hard-redirction to a specific
      function so that the current code work. Remove these
      hard-redirection after refactoring of the lower-level functions
      are done.

    - The argument “type” and “method” are only for the
      `log_pred_hurdle_brms` functions, not for Zresidual. The Zresidual
      method essential only repeat the calculation of multiple
      Z-residuals by using repeated random draws of U_i. But these
      arguments can be passed down to lower-level functions if it is
      called with Zresidual(model, …). The `type` should be renamed as
      `zr.type.hurdle`, and the `method` should be renamed as
      `sum.mcmc.method` which could be `post` or `iscv`.

    - The MPP method is coded as an option in the highest level
      Zresidual method with an argument called `randomized`. If
      `randomized=FALSE`, $U_i=0.5$, and `nrep` argument takes no effect
      as there is no need to replicate Z-residual.

- Organize the function list

  - Add @concept tags to each function, then the reference list will be
    organized. Example:

        #' @concept Zresidual
        #' @concept visualization
        #' @concept log_pred
        #' @concept test

  - add such lines to the \_pkgdown.yaml:

        reference:


          - title: "Visualization"
            desc: "Functions for visualizing entropy contributions and model noise."
            contents:
              - has_concept("visualization")

          - title: "Tabulation"
            desc: "Functions for formatting, printing, and exporting results tables."
            contents:
              - has_concept("tabulation")

    The pkgdown will then organize the functions by concept. A function
    can have multiple concepts.

- Refactor Bayesian workflow

  The function post_logrpp, iscv_logrpp etc need to be changed with U_i
  removed and merged into log_pred.brms, which will return only
  log_cdf/log_pmf (post or iscv), which will then be passed to Zresidual
  to generate replicated Z-residual with random draws of U_i. Once
  Zresidual method is modified to take only log_cdf and log_pmf as
  inputs, the following calculation is fast.

- Make generic function to compute Z-residual-based PPC p-values for
  Bayesian log_pred_mcmc object, which is an object of log_cdf and
  log_pmf for each MCMC samples
