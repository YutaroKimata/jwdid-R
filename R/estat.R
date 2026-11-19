#' Aggregate treatment effects
#'
#' Averages log ratios of counterfactual conditional means from a [jwdid()]
#' fit over observations and computes clustered delta-method standard errors.
#'
#' @param object A [jwdid()] fit.
#' @param type Aggregate post-baseline effects overall (`simple`), by cohort
#'   (`group`) or by period (`calendar`); aggregate by event time (`event`) or
#'   cohort-period (`attgt`); or pool all treated-cohort observations (`any`).
#'   Event aggregation includes pre-treatment effects when `never = TRUE`.
#' @param weights Nonnegative aggregation weights, as a vector aligned with the
#'   estimation data, a column name or a one-sided formula. Rows with zero
#'   weight are left out. Defaults to the estimation weights.
#' @param subset Logical vector, column name or one-sided formula selecting
#'   rows of the estimation data.
#' @param over Vector, column name or one-sided formula defining subgroups.
#' @param window Event-time range `c(lower, upper)` to keep, in time units.
#' @param cwindow Event-time range `c(lower, upper)` whose tails are collapsed
#'   into the endpoints. Use either `window` or `cwindow`.
#' @param asis Use the observed treatment indicator, instead of one, for the
#'   treated prediction.
#' @param level Confidence level.
#' @param pretrend For event aggregation with `never = TRUE`, test jointly that
#'   effects before the baseline are zero.
#'
#' @details
#' Each effect compares linear predictions with the treatment indicator at zero
#' and at one, holding covariates at their values. Group and calendar
#' aggregation leave out units with no period before their cohort in the
#' estimation sample. The pretrend test covers event times before the baseline.
#'
#' @return A data frame of class `jwdid_estat` with the grouping columns,
#'   `estimate`, `std.error`, `statistic`, `p.value`, `conf.low`, `conf.high`,
#'   `n` and `weight_sum`. Attributes hold the covariance (`vcov`), the
#'   gradient (`gradient`) and, if requested, the joint test (`pretrend`).
#'
#' @seealso [jwdid()], [plot.jwdid_estat()]
#' @examples
#' fit <- jwdid(trade, ivar = exporter^importer, tvar = year, gvar = first_year,
#'   fevar = ~ exporter^year + importer^year, never = TRUE, data = gravity)
#' es <- estat(fit, "event", pretrend = TRUE)
#' attr(es, "pretrend")
#' estat(fit, "event", window = c(-2, 3))
#' @export
estat <- function(object, type = c("simple", "group", "calendar", "event", "attgt", "any"),
                  weights = NULL, subset = NULL, over = NULL,
                  window = NULL, cwindow = NULL, asis = FALSE, level = 0.95, pretrend = FALSE) {
  type <- match.arg(type)
  if (!inherits(object, "jwdid")) stop("object must be a jwdid fit.")
  if (length(level) != 1L || !(level > 0 && level < 1)) stop("level must lie in (0, 1).")
  if (!is.null(window) && !is.null(cwindow)) stop("Use either window or cwindow.")
  bounds <- if (is.null(window)) cwindow else window
  if (!is.null(bounds) && (length(bounds) != 2L || bounds[1L] > bounds[2L])) {
    stop("A window is c(lower, upper) in time units.")
  }
  if (pretrend && (type != "event" || !object$never)) {
    stop("pretrend requires event aggregation and never = TRUE.")
  }

  d <- object$data
  w <- column(weights, d, d$.jw_w)
  if (any(!is.finite(w) | w < 0)) stop("Aggregation weights must be finite and nonnegative.")
  rows <- column(subset, d, rep(TRUE, nrow(d)))
  rows <- !is.na(rows) & rows & w > 0 & in_scope(d, type, object$never)
  event <- d$.jw_t - d$.jw_g
  if (!is.null(window)) rows <- rows & event >= window[1L] & event <= window[2L]
  if (!is.null(cwindow)) event <- pmin(pmax(event, cwindow[1L]), cwindow[2L])
  keys <- switch(type,
    simple = data.frame(aggregate = rep("simple", nrow(d))),
    any = data.frame(aggregate = rep("any", nrow(d))),
    group = data.frame(cohort = d$.jw_g),
    calendar = data.frame(time = d$.jw_t),
    event = data.frame(event = event),
    attgt = data.frame(cohort = d$.jw_g, time = d$.jw_t)
  )
  if (!is.null(over)) {
    keys$over <- column(over, d, NULL)
    rows <- rows & !is.na(keys$over)
  }
  if (!any(rows)) stop("No observations in the requested aggregation.")

  d <- d[rows, , drop = FALSE]
  w <- w[rows]
  keys <- keys[rows, , drop = FALSE]
  group <- do.call(interaction, c(keys, list(drop = TRUE, lex.order = TRUE)))
  weight_sum <- as.vector(rowsum(w, group))
  average <- Matrix::sparseMatrix(
    i = as.integer(group), j = seq_along(group), x = w / weight_sum[group],
    dims = c(nlevels(group), length(group))
  )

  b <- object$coefficients
  treatment <- jw_columns(d, names(b), object$design, "treatment")
  gradient <- as.matrix(average %*% ((if (asis) d$.jw_tr else 1) * treatment))
  dimnames(gradient) <- list(levels(group), names(b))
  estimate <- drop(gradient %*% b)
  covariance <- gradient %*% object$vcov %*% t(gradient)
  se <- sqrt(pmax(diag(covariance), 0))
  critical <- stats::qnorm((1 + level) / 2)

  out <- keys[match(seq_len(nlevels(group)), as.integer(group)), , drop = FALSE]
  rownames(out) <- NULL
  out$estimate <- estimate
  out$std.error <- se
  out$statistic <- estimate / se
  out$p.value <- 2 * stats::pnorm(-abs(out$statistic))
  out$conf.low <- estimate - critical * se
  out$conf.high <- estimate + critical * se
  out$n <- tabulate(group, nlevels(group))
  out$weight_sum <- weight_sum
  attr(out, "vcov") <- covariance
  attr(out, "gradient") <- gradient
  attr(out, "type") <- type
  if (pretrend) attr(out, "pretrend") <- pretrend_test(out, covariance, -object$gap * object$anticipation)
  class(out) <- c("jwdid_estat", "data.frame")
  out
}

# Rows each aggregation averages over.
in_scope <- function(d, type, never) {
  treated <- d$.jw_g > 0
  switch(type,
    simple = d$.jw_etr,
    group = ,
    calendar = d$.jw_etr & d$.jw_first < d$.jw_g,
    event = if (never) treated & d$.jw_tr != 0 else d$.jw_etr,
    attgt = ,
    any = treated
  )
}

pretrend_test <- function(out, covariance, baseline) {
  before <- which(out$event < baseline)
  if (!length(before)) stop("No pre-treatment effects in the requested window.")
  eig <- eigen(covariance[before, before, drop = FALSE], symmetric = TRUE)
  retain <- eig$values > max(eig$values, 0) * 1e-9
  rank <- sum(retain)
  if (!rank) stop("No estimable pre-treatment effects for the joint test.")
  chi2 <- sum(drop(crossprod(eig$vectors[, retain, drop = FALSE], out$estimate[before]))^2 /
    eig$values[retain])
  list(statistic = chi2, df = rank, p.value = stats::pchisq(chi2, rank, lower.tail = FALSE))
}

#' @export
print.jwdid_estat <- function(x, ...) {
  shown <- setdiff(names(x), c("statistic", "conf.low", "conf.high", "n", "weight_sum"))
  print.data.frame(x[shown], ...)
  invisible(x)
}

#' Plot aggregated treatment effects
#'
#' Plots event, calendar or cohort effects with confidence intervals using
#' [fixest::coefplot()].
#'
#' @param x A `jwdid_estat` table from event, calendar or group aggregation,
#'   with at most one `over` group.
#' @param xlab Horizontal axis label; defaults to the grouping column.
#' @param ylab Vertical axis label.
#' @param main Plot title.
#' @param ... Further arguments to [fixest::coefplot()].
#' @return `x`, invisibly.
#' @examples
#' fit <- jwdid(trade, ivar = exporter^importer, tvar = year, gvar = first_year,
#'   fevar = ~ exporter^year + importer^year, never = TRUE, data = gravity)
#' plot(estat(fit, "event"))
#' @importFrom graphics plot
#' @export
plot.jwdid_estat <- function(x, xlab = NULL, ylab = "Estimated treatment effect", main = "", ...) {
  key <- intersect(c("event", "time", "cohort"), names(x))[1L]
  if (is.na(key)) stop("Plots need event, calendar or group aggregation.")
  if ("over" %in% names(x)) stop("Subset to one over group before plotting.")
  ord <- order(x[[key]])
  fixest::coefplot(x$estimate[ord],
    x = x[[key]][ord], ci_low = x$conf.low[ord], ci_high = x$conf.high[ord],
    xlab = if (is.null(xlab)) key else xlab, ylab = ylab, main = main, ...
  )
  invisible(x)
}
