test_that("sparse columns reproduce fixest's design matrix", {
  d <- within(panel, {
    y <- y + 1
    t <- t / 10
    g <- g / 10
    constant <- id / 7
  })
  for (hettype in c("timecohort", "time", "cohort", "event", "twfe")) {
    for (cvar in list(NULL, ~x, ~ x + constant + I(2 * x))) {
      fit <- jwdid(y, id, t, g, cvar = cvar, hettype = hettype, data = d)
      probe <- fit$data
      for (tr in list(probe$.jw_tr, rep(1, nrow(probe)))) {
        probe$.jw_tr <- tr
        expected <- stats::model.matrix(fit$fixest, data = probe)
        actual <- as.matrix(jw_design_matrix(probe, colnames(expected), fit$design))
        expect_equal(actual, expected, tolerance = 1e-12, ignore_attr = TRUE, info = hettype)
      }
    }
  }
})
