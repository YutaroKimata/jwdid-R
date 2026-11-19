# Unit effects and their gradients from fixest's dense design matrix.
dense_effects <- function(fit) {
  b <- coef(fit)
  gradient <- stats::model.matrix(fit$fixest, data = fit$data)[, names(b)]
  gradient[, !startsWith(names(b), ".jw_cell::")] <- 0
  list(effect = drop(gradient %*% b), gradient = gradient)
}
