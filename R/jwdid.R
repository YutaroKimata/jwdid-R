#' Staggered difference-in-differences with PPML
#'
#' Estimates extended two-way fixed effects difference-in-differences models
#' by Poisson pseudo-maximum likelihood. Treatment cells are built from each
#' unit's first treatment period and estimated with [fixest::fepois()].
#'
#' @param yvar Nonnegative outcome, as a column name, string or expression.
#' @param ivar Unit identifier, as a column name or string. Combined
#'   identifiers can use `exporter^importer`.
#' @param tvar Numeric time variable, as a column name or string.
#' @param gvar First treatment period, as a column name or string. Zero marks
#'   never-treated units; the value must be constant within units.
#' @param cvar One-sided formula of covariates, such as `~ x1 + x2`.
#' @param fevar One-sided formula of additional fixed effects, such as
#'   `~ exporter^year + importer^year`.
#' @param never Use never-treated units as controls and estimate pre-treatment
#'   effects. Otherwise not-yet-treated observations serve as controls.
#' @param anticipation Positive integer: the baseline lies this many period
#'   steps before treatment.
#' @param hettype Treatment-effect heterogeneity: a separate effect for each
#'   cohort-period cell (`timecohort`), by period and pre/post status (`time`),
#'   by cohort and pre/post status (`cohort`), by time relative to treatment
#'   (`event`), or a single pre and post effect (`twfe`).
#' @param data A data frame.
#' @param subset Logical vector, column name or one-sided formula selecting
#'   rows of `data`.
#' @param weights Nonnegative probability weights, as a vector, column name or
#'   one-sided formula. Rows with zero weight are dropped.
#' @param offset Offset in the log mean, as a vector, column name or one-sided
#'   formula.
#' @param cluster One-sided clustering formula, such as
#'   `~ exporter + importer`. Defaults to the unit identifier.
#' @param ... Further arguments to [fixest::fepois()], such as `glm.tol`.
#'
#' @details
#' Units already treated in their first observed period are dropped. Without
#' never-treated units, the last cohort serves as the control group and
#' periods from its treatment onward are dropped. The period step is the
#' smallest spacing among the positive cohort and time values; a cohort first
#' treated at `g` has its baseline at `g - anticipation * step`.
#'
#' Rows with missing outcomes, covariates, identifiers or weights are dropped
#' before these rules, and rows with missing fixed effects, clusters or
#' offsets after them. Observations whose conditional mean is not identified
#' (separation) are detected and dropped before estimation.
#'
#' @return An object of class `jwdid`, with methods for [coef()], [vcov()],
#'   [nobs()], [fitted()], [summary()] and [estat()]. The underlying
#'   [fixest::fepois()] fit is the element `fixest`.
#'
#' @seealso [estat()]
#' @examples
#' fit <- jwdid(trade, ivar = exporter^importer, tvar = year, gvar = first_year,
#'   fevar = ~ exporter^year + importer^year, never = TRUE, data = gravity)
#' fit
#' estat(fit)
#' @export
jwdid <- function(yvar, ivar, tvar, gvar, cvar = NULL, fevar = NULL,
                  never = FALSE, anticipation = 1L,
                  hettype = c("timecohort", "time", "cohort", "event", "twfe"),
                  data, subset = NULL, weights = NULL, offset = NULL,
                  cluster = NULL, ...) {
  hettype <- match.arg(hettype)
  if (!is.numeric(anticipation) || length(anticipation) != 1L || is.na(anticipation) ||
    anticipation < 1 || anticipation != round(anticipation)) {
    stop("anticipation must be a positive integer.")
  }
  roles <- list(
    y = as_expression(substitute(yvar)), id = as_expression(substitute(ivar)),
    t = as_expression(substitute(tvar)), g = as_expression(substitute(gvar))
  )
  d <- jw_panel(data, roles, cvar, fevar, cluster, subset, weights, offset, parent.frame())

  timing <- jw_timing(d$.jw_id, d$.jw_t, d$.jw_g, anticipation)
  d <- d[timing$keep, , drop = FALSE]
  if (any(d$.jw_y < 0)) stop("PPML requires nonnegative outcomes.")
  cells <- jw_cells(d$.jw_t, d$.jw_g, timing$gap, timing$gmax, anticipation, never, hettype)
  d[names(cells)] <- cells
  if (nlevels(d$.jw_cell) < 2L) stop("No identified treatment cells.")
  d <- jw_covariates(d, cvar, hettype, timing$gap)
  design <- jw_design(d, timing$gmax)

  later <- grep("^\\.jw_(fe|cl)[0-9]+$|^\\.jw_off$", names(d), value = TRUE)
  d <- d[stats::complete.cases(d[later]) & is.finite(d$.jw_y) & is.finite(d$.jw_off), , drop = FALSE]
  separated <- jw_separation(d, design)
  dropped <- sort(d$.jw_row[separated])
  d <- d[!separated, , drop = FALSE]

  fit <- jw_fit(d, design, weighted = !is.null(weights), offset = !is.null(offset), ...)
  d <- d[fixest::obs(fit), , drop = FALSE]
  d$.jw_mu <- stats::fitted(fit)
  d$.jw_first <- stats::ave(d$.jw_t, d$.jw_id, FUN = min)

  structure(
    list(
      coefficients = stats::coef(fit), vcov = jw_vcov(fit, d, design), fixest = fit,
      data = d, design = design, never = never, hettype = hettype,
      gap = timing$gap, anticipation = anticipation, separated = dropped,
      call = match.call()
    ),
    class = "jwdid"
  )
}

jw_fit <- function(d, design, weighted, offset, ...) {
  clusters <- grep("^\\.jw_cl[0-9]+$", names(d), value = TRUE)
  cluster <- stats::reformulate(clusters)
  fit <- do.call(fixest::fepois, c(
    list(
      fml = design$formula, data = d,
      weights = if (weighted) ~.jw_w, offset = if (offset) ~.jw_off,
      vcov = if (length(clusters) == 1L) cluster else fixest::vcov_cluster(cluster, vcov_fix = FALSE),
      ssc = fixest::ssc(K.adj = FALSE, G.adj = TRUE, G.df = "min"),
      fixef.rm = "none"
    ),
    utils::modifyList(list(glm.tol = 1e-10, glm.iter = 10000L, fixef.tol = 1e-9), list(...))
  ))
  if (!isTRUE(fit$convStatus)) stop("PPML did not converge; no result returned.")
  fit
}
