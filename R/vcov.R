# Joint clustered covariance of the coefficients and the intercept left after
# absorbing fixed effects. Effects on the outcome scale depend on that
# intercept, so their delta-method variances need its covariances.

jw_vcov <- function(fit, d, design) {
  b <- stats::coef(fit)
  v <- stats::vcov(fit)
  x <- jw_design_matrix(d, names(b), design)
  q <- stats::fitted(fit) * d$.jw_w
  xm <- as.vector(Matrix::crossprod(x, q)) / sum(q)

  clusters <- lapply(d[grep("^\\.jw_cl[0-9]+$", names(d), value = TRUE)], fixest::to_integer)
  residual <- stats::residuals(fit, type = "response") * d$.jw_w / sum(q)
  scores <- fixest::estfun(fit)
  bread <- fixest::bread(fit) / stats::nobs(fit)
  cross <- numeric(length(b))
  variance <- 0
  # Inclusion-exclusion over cluster intersections, as in multiway clustering.
  for (k in seq_along(clusters)) {
    for (set in utils::combn(seq_along(clusters), k, simplify = FALSE)) {
      g <- do.call(fixest::to_integer, unname(clusters[set]))
      sign <- (-1)^(k + 1L)
      r <- drop(rowsum(residual, g, reorder = FALSE))
      cross <- cross + sign * drop(crossprod(rowsum(scores, g, reorder = FALSE), r))
      variance <- variance + sign * sum(r^2)
    }
  }
  size <- min(vapply(clusters, max, integer(1)))
  cross <- drop(crossprod(bread, cross)) * size / (size - 1L)
  vx <- drop(v %*% xm)
  covariance <- cross - vx
  joint <- rbind(
    cbind(v, `_cons` = covariance),
    `_cons` = c(covariance, variance * size / (size - 1L) - 2 * sum(xm * cross) + sum(xm * vx))
  )
  if (length(clusters) > 1L) joint <- clip_negative(joint, column_scales(x))
  joint
}

# Multiway covariances can be indefinite; clip negative eigenvalues of the
# covariance of standardized regressors.
clip_negative <- function(v, scales) {
  s <- c(pmax(scales, 1e-3), 1)
  eig <- eigen(v * outer(s, s), symmetric = TRUE)
  if (min(eig$values) >= -1e-10 * max(abs(eig$values))) return(v)
  warning("Multiway covariance was not positive semidefinite; eigenvalues clipped.")
  clipped <- eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors)) / outer(s, s)
  dimnames(clipped) <- dimnames(v)
  clipped
}

column_scales <- function(x) {
  n <- nrow(x)
  means <- Matrix::colMeans(x)
  sqrt(pmax(Matrix::colSums(x^2) - n * means^2, 0) / (n - 1L))
}
