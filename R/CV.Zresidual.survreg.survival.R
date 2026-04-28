#' @rdname CV.Zresidual
#' @method CV.Zresidual survreg
#' @export
CV.Zresidual.survreg <- function(object,
                                 nfolds = NULL,
                                 foldlist = NULL,
                                 data = NULL,
                                 nrep = 1,
                                 log_pointpred = NULL,
                                 mcmc_summarize = c("post", "iscv"),
                                 type = NULL,
                                 randomized = TRUE,
                                 seed = NULL,
                                 ...) {
  .cv_run_survival_cv(
    object = object,
    nfolds = nfolds,
    foldlist = foldlist,
    data = data,
    nrep = nrep,
    log_pointpred = log_pointpred,
    mcmc_summarize = mcmc_summarize,
    type = type,
    randomized = randomized,
    seed = seed,
    ...
  )
}
