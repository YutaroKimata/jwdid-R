# Multiway clustered covariances can be indefinite; negative eigenvalues of the
# covariance of standardized regressors are clipped.

jw_vcov <- function(fit, d, design) {
  v <- stats::vcov(fit)
  if (length(grep("^\\.jw_cl[0-9]+$", names(d))) < 2L) return(v)
  x <- jw_design_matrix(d, names(stats::coef(fit)), design)
  clip_negative(v, column_scales(x))
}

clip_negative <- function(v, scales) {
  s <- pmax(scales, 1e-3)
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
