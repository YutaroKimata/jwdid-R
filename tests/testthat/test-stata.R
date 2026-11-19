test_that("aggregates match the saved Stata results", {
  reference <- read.csv(test_path("fixtures", "stata.csv"))
  last_cohort <- within(panel, g[g == 0] <- 8)
  cases <- list(
    ppml_never = list(yvar = "y"),
    ppml_covariates = list(yvar = "y", cvar = ~x),
    ppml_notyet = list(yvar = "y", cvar = ~x, never = FALSE),
    ppml_no_never = list(yvar = "y", cvar = ~x, never = FALSE, data = last_cohort),
    gravity = list(
      yvar = "trade", ivar = "exporter^importer", tvar = "year",
      fevar = ~ exporter^year + importer^year, data = trade
    ),
    anticipation = list(yvar = "y", cvar = ~x, anticipation = 2L),
    window = list(yvar = "y", cvar = ~x),
    cwindow = list(yvar = "y", cvar = ~x)
  )
  checked <- 0L
  for (case in names(cases)) {
    args <- list(ivar = "id", tvar = "t", gvar = "g", never = TRUE, data = panel)
    args[names(cases[[case]])] <- cases[[case]]
    fit <- do.call(jwdid, args)
    expected <- reference[reference$case == case, ]
    for (type in unique(expected$type)) {
      result <- estat(fit, type,
        window = if (case == "window") c(-2, 3), cwindow = if (case == "cwindow") c(-2, 3))
      target <- expected[expected$type == type, ]
      expect_equal(result$estimate, target$Stata_estimate, tolerance = 1e-6, info = paste(case, type))
      expect_equal(result$std.error, target$Stata_se, tolerance = 1e-6, info = paste(case, type))
      checked <- checked + nrow(result)
    }
  }
  expect_equal(checked, nrow(reference))
})
