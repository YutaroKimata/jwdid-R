# Unit effects and their gradients from fixest's dense design matrices.
dense_effects <- function(fit) {
  e <- fit$data
  b <- coef(fit)
  design <- function(tr) stats::model.matrix(fit$fixest, data = transform(e, .jw_tr = tr))[, names(b)]
  gradient <- design(1) - design(0)
  list(effect = drop(gradient %*% b), gradient = gradient)
}
