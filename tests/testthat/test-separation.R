test_that("separated observations are detected and dropped", {
  patterns <- list(
    unit = panel$id %in% c(1, 2),
    post_cell = panel$g == 3 & panel$t == 5,
    base_period = panel$g == 3 & panel$t == 2
  )
  for (pattern in names(patterns)) {
    d <- panel
    d$y[patterns[[pattern]]] <- 0
    fit <- jwdid(y, id, t, g, never = TRUE, data = d)
    expect_equal(fit$separated, which(patterns[[pattern]]), info = pattern)
  }
})

test_that("covariates do not hide base-period separation", {
  d <- panel
  zero <- d$g == 3 & d$t == 2
  d$y[zero] <- 0
  fit <- jwdid(y, id, t, g, cvar = ~x, never = TRUE, data = d)
  expect_true(all(which(zero) %in% fit$separated))
  expect_lt(max(abs(coef(fit))), 5)
})
