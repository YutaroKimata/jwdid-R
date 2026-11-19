#' @exportS3Method stats::coef
coef.jwdid <- function(object, ...) object$coefficients

#' @exportS3Method stats::vcov
vcov.jwdid <- function(object, ...) object$vcov

#' @exportS3Method stats::nobs
nobs.jwdid <- function(object, ...) nrow(object$data)

#' @exportS3Method stats::fitted
fitted.jwdid <- function(object, ...) stats::fitted(object$fixest)

#' @export
print.jwdid <- function(x, ...) {
  cat("Extended TWFE difference-in-differences (PPML)\n")
  cat(sprintf(
    "Observations: %d   Clusters: %d   Treatment cells: %d   Separated: %d\n",
    stats::nobs(x), length(unique(x$data$.jw_cl1)), nlevels(x$data$.jw_cell) - 1L, length(x$separated)
  ))
  cat(sprintf("Controls: %s\n", if (x$never) "never treated" else "not yet treated"))
  cat("Use estat() for aggregated effects and summary() for coefficients.\n")
  invisible(x)
}

#' @export
summary.jwdid <- function(object, ...) {
  b <- object$coefficients
  se <- sqrt(diag(stats::vcov(object)))
  z <- b / se
  structure(
    list(coefficients = cbind(Estimate = b, `Std. Error` = se, `z value` = z,
      `Pr(>|z|)` = 2 * stats::pnorm(-abs(z))), fit = object),
    class = "summary.jwdid"
  )
}

#' @export
print.summary.jwdid <- function(x, ...) {
  print(x$fit)
  cat("\n")
  stats::printCoefmat(x$coefficients, ...)
  invisible(x)
}
