test_that("aggregates match dense counterfactual calculations", {
  d <- within(panel, {
    x2 <- sin(id / 2 + t)
    x3 <- cos(id / 7 + t / 3)
    w <- 1 + id %% 4 / 4
    off <- sin(t) / 8
  })
  fit <- jwdid(y, id, t, g, cvar = ~ x + x2 + x3, never = TRUE, weights = ~w, offset = ~off, data = d)
  e <- fit$data
  influence <- rowsum(fixest::estfun(fit$fixest), e$.jw_id) %*% (fixest::bread(fit$fixest) / nobs(fit))
  clusters <- nrow(influence)
  covariance <- crossprod(influence) * clusters / (clusters - 1)

  keep <- e$.jw_g > 0 & e$.jw_t >= 2
  event <- pmax(-2, pmin(3, e$.jw_t - e$.jw_g))
  over <- e$t %% 2
  weights <- e$w
  result <- estat(fit, "event", subset = ~ t >= 2, over = ~ t %% 2, cwindow = c(-2, 3))
  unit <- dense_effects(fit)
  groups <- lapply(seq_len(nrow(result)), function(i) {
    which(keep & event == result$event[i] & over == result$over[i])
  })
  expected <- t(vapply(groups, function(i) colSums(unit$gradient[i, , drop = FALSE] * weights[i]) / sum(weights[i]),
    numeric(ncol(unit$gradient))))
  expect_equal(unname(attr(result, "gradient")), unname(expected), tolerance = 1e-10)
  expect_equal(result$estimate, vapply(groups, function(i) stats::weighted.mean(unit$effect[i], weights[i]), 1),
    tolerance = 1e-10)
  expect_equal(unname(attr(result, "vcov")), unname(expected %*% covariance %*% t(expected)), tolerance = 1e-10)
})

test_that("group and calendar effects leave out units without a pre-treatment period", {
  fit <- jwdid(y, id, t, g, never = TRUE, data = within(panel, y[g == 3 & t <= 2] <- 0))
  e <- fit$data
  expect_false(any(e$.jw_g == 3 & e$.jw_t <= 2))
  expect_false(3 %in% estat(fit, "group")$cohort)
  expect_equal(sum(estat(fit, "calendar")$n), sum(e$.jw_etr & e$.jw_g != 3))
  expect_equal(estat(fit, "simple")$n, sum(e$.jw_etr))
})

test_that("event and overall effects recover the true log effects", {
  oracle <- jwdid(mu, ivar = exporter^importer, tvar = year, gvar = first_year,
    fevar = ~ exporter^year + importer^year, never = TRUE, data = gravity,
    glm.tol = 1e-8, glm.iter = 100)
  treated <- gravity[gravity$first_year > 0, ]
  truth <- tapply(treated$tau, treated$year - treated$first_year, mean)
  event <- estat(oracle, "event")
  expect_equal(event$estimate, as.vector(truth[as.character(event$event)]), tolerance = 1e-6)
  expect_equal(estat(oracle)$estimate, mean(gravity$tau[gravity$rta == 1]), tolerance = 1e-6)
})

test_that("the baseline is a zero effect and the pretrend test is defined", {
  fit <- jwdid(y, id, t, g, cvar = ~x, never = TRUE, data = panel)
  baseline <- estat(fit, "event", window = c(-1, -1))
  expect_true(all(baseline$estimate == 0 & baseline$std.error == 0))
  test <- attr(estat(fit, "event", pretrend = TRUE), "pretrend")
  expect_equal(test$df, 5)
  expect_true(test$p.value > 0 && test$p.value < 1)
  expect_error(estat(fit, "simple", pretrend = TRUE), "pretrend")
})

test_that("tables print compactly and plot", {
  fit <- jwdid(y, id, t, g, never = TRUE, data = panel)
  shown <- capture.output(print(estat(fit)))
  expect_false(any(grepl("weight_sum|statistic|conf.low", shown)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (type in c("event", "calendar", "group")) {
    result <- estat(fit, type)
    expect_identical(plot(result), result)
  }
})
