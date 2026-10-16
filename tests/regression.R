library(jwdid.r)
stopifnot(!"package:fixest" %in% search())
rd <- tools::Rd_db("jwdid.r")[["jwdid.Rd"]]
arguments <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\arguments"), rd)[[1L]]
items <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\item"), arguments)
documented <- vapply(items, function(x) paste(unlist(x[[1L]]), collapse = ""), character(1))
stopifnot(identical(names(formals(jwdid)), documented))
panel <- read.csv("fixtures/panel.csv")
gravity <- read.csv("fixtures/gravity.csv")
reference <- read.csv("fixtures/stata.csv")
reference$predict <- unname(c(link = "xb", response = "mu")[reference$predict])
stopifnot(!anyNA(reference$predict))
last_cohort <- panel
last_cohort$g[last_cohort$g == 0] <- 8
cases <- list(
  ppml_never = list(yvar = "y"),
  ppml_covariates = list(yvar = "y", cvar = ~x),
  ppml_notyet = list(yvar = "y", cvar = ~x, never = FALSE),
  ppml_no_never = list(yvar = "y", cvar = ~x, never = FALSE, data = last_cohort),
  gravity = list(yvar = "trade", ivar = "exporter^importer", tvar = "year",
    fevar = ~exporter^year + importer^year, data = gravity),
  anticipation = list(yvar = "y", cvar = ~x, anticipation = 2L),
  window = list(yvar = "y", cvar = ~x),
  cwindow = list(yvar = "y", cvar = ~x)
)
fits <- list()
checked <- 0L
for (case in names(cases)) {
  args <- list(ivar = "id", tvar = "t", gvar = "g", never = TRUE, data = panel)
  args[names(cases[[case]])] <- cases[[case]]
  fits[[case]] <- do.call(jwdid, args)
  expected <- reference[reference$case == case, ]
  for (type in unique(expected$type)) {
    for (scale in unique(expected$predict)) {
      options <- list(object = fits[[case]], type = type, predict = scale)
      if (case == "window") options$window <- c(-2, 3)
      if (case == "cwindow") options$cwindow <- c(-2, 3)
      result <- do.call(estat, options)
      target <- expected[expected$type == type & expected$predict == scale, ]
      stopifnot(
        nrow(result) == nrow(target),
        all(is.finite(result$estimate)), all(is.finite(result$std.error)),
        max(abs(result$estimate - target$Stata_estimate)) < 1e-5,
        max(abs(result$std.error - target$Stata_se)) < 1e-5
      )
      checked <- checked + nrow(result)
    }
  }
}
stopifnot(checked == nrow(reference))

explicit <- jwdid(yvar = trade, ivar = id, tvar = year, gvar = g,
  fevar = ~exporter^year + importer^year, never = TRUE, data = gravity)
bare <- jwdid(trade, exporter^importer, year, g, NULL,
  ~exporter^year + importer^year, never = TRUE, data = gravity)
default <- jwdid(y, id, t, g, cvar = ~x, data = panel)
combined <- fits$gravity
stopifnot(
  isTRUE(all.equal(coef(default), coef(fits$ppml_notyet), tolerance = 1e-12)),
  isTRUE(all.equal(coef(combined), coef(explicit), tolerance = 1e-12)),
  isTRUE(all.equal(coef(combined), coef(bare), tolerance = 1e-12)),
  isTRUE(all.equal(vcov(combined), vcov(bare), tolerance = 1e-12)),
  identical(estat(combined), estat(combined, predict = "xb")),
  inherits(estat(combined), "jwdid_estat"),
  identical(class(combined), "fixest"),
  is.null(combined$jwdid$stata),
  is.null(combined$model_info$offset), is.null(combined$model_info$weights),
  isTRUE(all.equal(predict(combined), fitted(combined))),
  isTRUE(all.equal(coef(summary(combined)), coef(combined))),
  isTRUE(all.equal(fitted(combined), fitted(explicit), tolerance = 1e-12)),
  isTRUE(all.equal(vcov(combined), vcov(explicit), tolerance = 1e-12)),
  identical(combined$jwdid$source_rows, explicit$jwdid$source_rows)
)
restored <- unserialize(serialize(combined, NULL))
stopifnot(identical(estat(restored, "event"), estat(combined, "event")),
  identical(estat(restored, "event", predict = "mu"), estat(combined, "event", predict = "mu")))
shown <- capture.output(print(estat(combined)))
stopifnot(!any(grepl("weight_sum|statistic|conf.low", shown)),
  identical(getS3method("print", "jwdid_estat"), jwdid.r:::print.jwdid_estat))
event <- estat(fits$ppml_covariates, "event", pretrend = TRUE)
stopifnot(is.finite(attr(event, "pretrend")$p.value))
stopifnot(identical(getS3method("plot", "jwdid_estat"), plot.jwdid_estat))
pdf(file = NULL)
for (type in c("event", "calendar", "group")) {
  result <- estat(combined, type)
  stopifnot(identical(plot(result, pch = 19), result))
}
dev.off()
cat(checked, "ATT estimates and standard errors match saved Stata results.\n")

example <- local({
  data(gravity)
  gravity
})
rate <- c(`2003` = 0.06, `2006` = 0.04, `2009` = 0.02, `0` = 0)
stopifnot(
  nrow(example) == 16L * 15L * 12L,
  any(example$trade == 0), all(example$trade >= 0), all(example$mu > 0),
  all(example$rta == (example$first_year > 0 & example$year >= example$first_year)),
  isTRUE(all.equal(example$tau, unname(rate[as.character(example$first_year)]) *
    pmax(example$year - example$first_year + 1, 0) * example$rta))
)
oracle <- jwdid(mu, ivar = exporter^importer, tvar = year, gvar = first_year,
  fevar = ~exporter^year + importer^year, never = TRUE, data = example,
  glm.tol = 1e-8, glm.iter = 100)
truth <- aggregate(tau ~ event,
  transform(example[example$first_year > 0, ], event = year - first_year), mean)
event <- estat(oracle, "event")
stopifnot(
  max(abs(event$estimate - truth$tau[match(event$event, truth$event)])) < 1e-6,
  abs(estat(oracle)$estimate - mean(example$tau[example$rta == 1])) < 1e-6
)
cat("PPML DGP: event and overall log effects match their known values.\n")

local({
  panel <- read.csv("fixtures/panel.csv")
  for (zero in list(panel$id %in% c(1, 2), panel$g == 3 & panel$t == 5)) {
    d <- panel
    d$y[zero] <- 0
    fit <- jwdid(y, id, t, g, never = TRUE, data = d)
    stopifnot(identical(fit$jwdid$separated_rows, which(zero)))
  }
  panel$x2 <- sin(panel$id / 2 + panel$t)
  panel$x3 <- cos(panel$id / 7 + panel$t / 3)
  panel$w <- 1 + panel$id %% 4 / 4
  panel$off <- sin(panel$t) / 8
  fit <- jwdid(y, id, t, g, cvar = ~x + x2 + x3,
    never = TRUE, weights = ~w, offset = ~off, data = panel)
  stopifnot(identical(parent.env(fit$call_env), environment()))
  d <- fit$jwdid$data
  original_treatment <- d$.jw_tr
  keep <- d$.jw_g > 0 & d$.jw_tr != 0 & d$.jw_t >= 2
  event <- pmax(-2, pmin(3, d$.jw_event))
  over <- d$.jw_t %% 2
  weights <- 1 + abs(d$.jw_raw001)
  d0 <- d1 <- d
  d0$.jw_tr <- 0
  d1$.jw_tr <- 1
  x0 <- stats::model.matrix(fit, data = d0)
  x1 <- stats::model.matrix(fit, data = d1)
  x <- stats::model.matrix(fit, data = d)
  b <- coef(fit)
  q <- fitted(fit) * fit$jwdid$weights
  xm <- drop(crossprod(x, q)) / sum(q)
  units <- fit$jwdid$data$id
  influence_b <- rowsum(fixest::estfun(fit), units) %*% (fixest::bread(fit) / nobs(fit))
  influence_c <- drop(rowsum(residuals(fit, type = "response") * fit$jwdid$weights, units)) / sum(q) -
    drop(influence_b %*% xm)
  clusters <- length(unique(units))
  joint_vcov <- crossprod(cbind(influence_b, influence_c)) * clusters / (clusters - 1L)

  for (scale in c("xb", "mu")) {
    result <- estat(fit, "event", predict = scale, weights = weights,
      subset = ~.jw_t >= 2, over = over, cwindow = c(-2, 3))
    if (scale == "xb") {
      gradient <- cbind(x1 - x0, `_cons` = 0)
      effects <- drop((x1 - x0) %*% b)
    } else {
      eta <- log(fitted(fit))
      p0 <- exp(eta + drop((x0 - x) %*% b))
      p1 <- exp(eta + drop((x1 - x) %*% b))
      effects <- p1 - p0
      gradient <- cbind(x1 * p1 - x0 * p0, `_cons` = effects)
    }
    groups <- lapply(seq_len(nrow(result)), function(i) {
      which(keep & event == result$event[i] & over == result$over[i])
    })
    expected_gradient <- t(vapply(groups, function(idx) {
      colSums(gradient[idx, , drop = FALSE] * weights[idx]) / sum(weights[idx])
    }, numeric(length(b) + 1L)))
    expected_effects <- vapply(groups, function(idx) {
      weighted.mean(effects[idx], weights[idx])
    }, numeric(1))
    expected_vcov <- expected_gradient %*% joint_vcov %*% t(expected_gradient)
    stopifnot(
      isTRUE(all.equal(unname(attr(result, "gradient")), unname(expected_gradient), tolerance = 1e-10)),
      isTRUE(all.equal(unname(result$estimate), expected_effects, tolerance = 1e-10)),
      isTRUE(all.equal(unname(attr(result, "vcov")), unname(expected_vcov), tolerance = 1e-10)),
      isTRUE(all.equal(unname(result$std.error), sqrt(pmax(diag(expected_vcov), 0)), tolerance = 1e-10)),
      identical(fit$jwdid$data$.jw_tr, original_treatment)
    )
    baseline <- estat(fit, "event", predict = scale, window = c(-1, -1))
    stopifnot(all(baseline$estimate == 0), all(baseline$std.error == 0),
      all(attr(baseline, "gradient") == 0))
  }
  cat("ATT gradients match dense counterfactual calculations.\n")
})

local({
  d <- panel
  d$y <- d$y + 1
  d$t <- d$t / 10
  d$g <- d$g / 10
  d$constant <- d$id / 7
  for (hettype in c("timecohort", "time", "cohort", "event", "twfe")) {
    for (cvar in list(NULL, ~1, ~x + constant + I(2 * x))) {
      fit <- jwdid(y, id, t, g, cvar = cvar, hettype = hettype, data = d)
      probe <- fit$jwdid$data
      design <- jwdid.r:::.jwdid_design(probe, names(coef(fit)))
      for (treatment in list(probe$.jw_tr, rep(1, nrow(probe)))) {
        probe$.jw_tr <- treatment
        expected <- stats::model.matrix(fit, data = probe)
        stopifnot(identical(colnames(expected), names(coef(fit))))
        actual <- matrix(0, nrow(expected), ncol(expected), dimnames = dimnames(expected))
        for (block in design) {
          actual[cbind(block$row, block$column)] <- block$value *
            if (block$treatment) treatment[block$row] else 1
        }
        stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-12, check.attributes = FALSE)))
      }
    }
  }
  cat("Sparse design matches model.matrix for every heterogeneity specification.\n")
})
