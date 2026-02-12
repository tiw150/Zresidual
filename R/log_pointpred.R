log_pointpred <- function(fit, data = NULL, config = list(), ...) {
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  
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
      m <- tryCatch(getS3method(generic, cls1, optional = TRUE), error = function(e) NULL)
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
  
  config <- config %||% list()
  type <- lp_norm_key(config$type %||% "")
  
  backend <- lp_norm_key(config$backend %||% lp_backend_from_fit(fit))
  family  <- lp_norm_key(config$family  %||% lp_family_from_fit(fit))
  
  base <- paste0("log_pointpred_", backend, "_", family)
  
  fname <- if (type == "") {
    base
  } else if (type == "count") {
    paste0(base, "_count")
  } else if (type == "zero") {
    if (!grepl("^hurdle(_|$)", family)) {
      stop("config$type='zero' is only valid for hurdle families.", call. = FALSE)
    }
    paste0("log_pointpred_", backend, "_hurdle_zero")
  } else {
    paste0(base, "_", type)
  }
  
  ns <- environment(log_pointpred)
  f <- get0(fname, envir = ns, inherits = TRUE)
  
  if (is.null(f)) {
    stop(sprintf("Cannot find underlying function '%s'.", fname), call. = FALSE)
  }
  
  f(fit, data = data, ...)
}
