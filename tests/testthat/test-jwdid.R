test_that("identifiers can be bare names, strings or combined", {
  bare <- jwdid(trade, exporter^importer, year, g, NULL, ~ exporter^year + importer^year,
    never = TRUE, data = trade)
  quoted <- jwdid("trade", "exporter^importer", "year", "g",
    fevar = ~ exporter^year + importer^year, never = TRUE, data = trade)
  single <- jwdid(trade, id, year, g, fevar = ~ exporter^year + importer^year, never = TRUE, data = trade)
  expect_equal(coef(bare), coef(quoted))
  expect_equal(coef(bare), coef(single))
  expect_equal(vcov(bare), vcov(single))
})

test_that("the not-yet-treated design is the default", {
  expect_equal(
    coef(jwdid(y, id, t, g, cvar = ~x, data = panel)),
    coef(jwdid(y, id, t, g, cvar = ~x, never = FALSE, data = panel))
  )
})

test_that("zero weights drop rows", {
  d <- within(panel, w <- ifelse(id %% 10 == 0, 0, 1 + id %% 4 / 4))
  expect_equal(
    coef(jwdid(y, id, t, g, never = TRUE, weights = ~w, data = d)),
    coef(jwdid(y, id, t, g, never = TRUE, weights = ~w, data = d[d$w > 0, ]))
  )
})

test_that("missing fixed effects drop rows after the sample rules", {
  units <- unique(panel$id[panel$g == 3])[1:5]
  d <- within(panel, {
    fe <- id %% 5
    fe[id %in% units & t <= 2] <- NA
  })
  fit <- jwdid(y, id, t, g, fevar = ~fe, never = TRUE, data = d)
  expect_equal(nobs(fit), nrow(panel) - 10)
})

test_that("methods summarize the fit", {
  fit <- jwdid(y, id, t, g, never = TRUE, data = panel)
  expect_s3_class(fit, "jwdid")
  expect_equal(nobs(fit), nrow(panel))
  expect_equal(fitted(fit), unname(stats::fitted(fit$fixest)))
  expect_equal(unname(summary(fit)$coefficients[, 2]), unname(sqrt(diag(vcov(fit)))))
  expect_output(print(fit), "Observations: 1440")
})

test_that("invalid inputs stop with a reason", {
  expect_error(jwdid(y, id, t, g, data = transform(panel, .jw_a = 1)), "reserved")
  expect_error(jwdid(y, id, t, g, data = rbind(panel, panel[1, ])), "Duplicate")
  expect_error(jwdid(y, id, t, g, data = transform(panel, y = -y)), "nonnegative")
  expect_error(jwdid(y, id, t, g, data = transform(panel, g = ifelse(t == 1, 0, g))), "constant")
  expect_error(jwdid(y, id, t, g, weights = ~ -1 + 0 * t, data = panel), "nonnegative")
  expect_error(jwdid(y, id, t, g, anticipation = 0, data = panel), "anticipation")
})
