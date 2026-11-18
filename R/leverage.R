# Residuals of each unit scaled by the inverse of one minus their leverage on
# the regressors and the non-unit fixed effects.

leverage_residuals <- function(fit, d, x) {
  lambda <- stats::fitted(fit) * d$.jw_w
  e <- stats::residuals(fit, type = "response")
  unit <- fixest::to_integer(d$.jw_id)
  others <- c(".jw_t", grep("^\\.jw_fe[0-9]+$", names(d), value = TRUE))
  groups <- lapply(d[c(".jw_id", others)], fixest::to_integer)

  xt <- fixest::demean(as.matrix(x), f = groups, weights = lambda, tol = 1e-12, iter = 10000L,
    notes = FALSE)
  xv <- xt %*% solve(crossprod(xt, xt * lambda))

  dummies <- do.call(cbind, lapply(groups[-1L], function(g) {
    Matrix::sparseMatrix(i = seq_along(g), j = g, x = 1, dims = c(length(g), max(g)))
  }))
  units <- Matrix::sparseMatrix(i = seq_along(unit), j = unit, x = 1)
  theta <- lambda / rowsum(lambda, unit)[unit]
  dt <- dummies - units %*% Matrix::crossprod(units, dummies * theta)
  hessian <- eigen(as.matrix(Matrix::crossprod(dt, dt * lambda)), symmetric = TRUE)
  keep <- hessian$values > max(hessian$values) * 1e-10
  vphi <- hessian$vectors[, keep] %*% (t(hessian$vectors[, keep]) / hessian$values[keep])

  entries <- Matrix::summary(methods::as(dt, "TsparseMatrix"))
  entries <- split(entries[c("i", "j", "x")], unit[entries$i])
  rows <- split(seq_along(unit), unit)
  adjusted <- numeric(length(e))
  for (g in seq_along(rows)) {
    r <- rows[[g]]
    cell <- entries[[as.character(g)]]
    cols <- unique(cell$j)
    dg <- matrix(0, length(r), length(cols))
    dg[cbind(match(cell$i, r), match(cell$j, cols))] <- cell$x
    lev <- lambda[r] * (xv[r, , drop = FALSE] %*% t(xt[r, , drop = FALSE]) +
      dg %*% vphi[cols, cols, drop = FALSE] %*% t(dg))
    adjusted[r] <- pseudo_solve(diag(length(r)) - lev, e[r])
  }
  list(demeaned = xt, residual = adjusted)
}

# Units with leverage one in some direction have zero residuals there.
pseudo_solve <- function(a, b) {
  s <- svd(a)
  keep <- s$d > max(s$d) * 1e-10
  drop(s$v[, keep, drop = FALSE] %*% (crossprod(s$u[, keep, drop = FALSE], b) / s$d[keep]))
}
