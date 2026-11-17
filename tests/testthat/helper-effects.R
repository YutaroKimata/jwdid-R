# Unit effects and their gradients from fixest's dense design matrices.
dense_effects <- function(fit, scale) {
  e <- fit$data
  b <- coef(fit)
  design <- function(tr) stats::model.matrix(fit$fixest, data = transform(e, .jw_tr = tr))[, names(b)]
  x <- design(e$.jw_tr)
  x0 <- design(0)
  x1 <- design(1)
  if (scale == "xb") return(list(effect = drop((x1 - x0) %*% b), gradient = cbind(x1 - x0, 0)))
  p0 <- fitted(fit) * exp(drop((x0 - x) %*% b))
  p1 <- fitted(fit) * exp(drop((x1 - x) %*% b))
  list(effect = p1 - p0, gradient = cbind(x1 * p1 - x0 * p0, p1 - p0))
}
