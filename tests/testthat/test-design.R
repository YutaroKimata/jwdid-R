test_that("sparse columns reproduce fixest's design matrix", {
  d <- within(panel, {
    y <- y + 1
    t <- t / 10
    g <- g / 10
    constant <- id / 7
  })
  for (cvar in list(NULL, ~x, ~ x + constant + I(2 * x))) {
    fit <- jwdid(y, id, t, g, cvar = cvar, data = d)
    expected <- stats::model.matrix(fit$fixest, data = fit$data)
    actual <- as.matrix(jw_design_matrix(fit$data, colnames(expected), fit$design))
    expect_equal(actual, expected, tolerance = 1e-12, ignore_attr = TRUE)
  }
})
