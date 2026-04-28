
log_pointpred <- function(fit, data, type = NULL, ...) {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  
  if (missing(data) || is.null(data)) {
    stop("log_pointpred: `data` must be provided.", call. = FALSE)
  }
  
  lp_norm_key <- function(x) {
    x <- tolower(as.character(x %||% ""))
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }
  
  lp_backend_from_fit <- function(fit) {
    cls <- class(fit)
    
    pkg <- attr(cls, "package")
    if (!is.null(pkg) && length(pkg) && nzchar(pkg[1])) return(lp_norm_key(pkg[1]))
    
    cls1 <- cls[1]
    for (generic in c("predict", "fitted", "logLik", "print")) {
      m <- tryCatch(utils::getS3method(generic, cls1, optional = TRUE), error = function(e) NULL)
      if (is.function(m)) {
        ns <- environmentName(environment(m))
        if (nzchar(ns) && ns != "R_GlobalEnv") return(lp_norm_key(ns))
      }
    }
    
    lp_norm_key(sub("fit$", "", cls1, ignore.case = TRUE))
  }
  
  lp_family_from_fit <- function(fit) {
    fam <- NULL
    if (!is.null(fit$family)) fam <- fit$family$family %||% fit$family
    if (is.null(fam) && !is.null(fit$distribution)) fam <- fit$distribution
    if (is.null(fam) && !is.null(attr(fit, "family"))) fam <- attr(fit, "family")
    if (is.null(fam)) fam <- class(fit)[1]
    lp_norm_key(fam)
  }
  
  lp_type_from_fit <- function(fit) {
    if (inherits(fit, "coxph")) {
      frailty_terms <- tryCatch(attr(fit$terms, "specials")$frailty, error = function(e) NULL)
      if (!is.null(frailty_terms) && length(frailty_terms) > 0L) {
        return("frailty")
      }
    }
    NULL
  }
  
  backend <- lp_backend_from_fit(fit)
  family  <- lp_family_from_fit(fit)
  type_key <- lp_norm_key(type %||% lp_type_from_fit(fit) %||% "")
  
  base <- paste0("log_pointpred_", backend, "_", family)
  
  fname <- if (type_key == "") {
    base
  } else if (type_key == "count") {
    paste0(base, "_count")
  } else if (type_key == "zero") {
    if (!grepl("^hurdle(_|$)", family)) {
      stop("type='zero' is only valid for hurdle families.", call. = FALSE)
    }
    paste0("log_pointpred_", backend, "_hurdle_zero")
  } else {
    paste0(base, "_", type_key)
  }
  
  f <- get0(fname, envir = .GlobalEnv, inherits = FALSE)
  if (is.null(f)) {
    f <- get0(fname, envir = environment(log_pointpred), inherits = TRUE)
  }
  
  if (is.null(f) || !is.function(f)) {
    stop(
      sprintf(
        paste0(
          "Cannot find underlying function '%s' in .GlobalEnv or package namespace.\n",
          "Hint: implement '%s' using the pattern log_pointpred_<backend>_<family>[_<type>]."
        ),
        fname, fname
      ),
      call. = FALSE
    )
  }
  
  f(fit, data = data, ...)
}
