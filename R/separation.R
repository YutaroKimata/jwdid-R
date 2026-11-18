# Separation: observations with y = 0 whose conditional mean the model can
# drive to zero, so the PPML estimate does not exist. Each pass drops zero-sum
# fixed-effect groups, then runs the ReLU check of Correia, Guimaraes and
# Zylkin; passes repeat until nothing more separates.

jw_separation <- function(d, design) {
  separated <- logical(nrow(d))
  fixed <- c(".jw_id", ".jw_t", grep("^\\.jw_fe[0-9]+$", names(d), value = TRUE))
  repeat {
    live <- d[!separated, , drop = FALSE]
    groups <- lapply(live[fixed], fixest::to_integer)
    y <- live$.jw_y
    found <- Reduce(`|`, lapply(groups, function(g) rowsum(y, g)[g, 1L] == 0))
    if (!any(found) && any(y == 0)) found <- relu_separation(live, design, groups)
    if (!any(found)) return(separated)
    if (all(found)) stop("No observations remain after separation removal.")
    separated[which(!separated)[found]] <- TRUE
  }
}

relu_separation <- function(d, design, groups) {
  boundary <- d$.jw_y == 0
  # Cells stay regressors: absorbing them as fixed effects stalls the projection.
  x <- suppressWarnings(as.matrix(jw_design_matrix(d, design$columns, design)))
  relu_iterate(boundary, projection(x, groups, boundary))
}

# Least-squares residual of `u` on the fixed effects and regressors, with the
# fit held at `u` where y > 0 by feeding the violation back into the target.
projection <- function(x, groups, boundary, weight = 100, maxiter = 1000L) {
  for (j in seq_len(ncol(x))) {
    scale <- sqrt(sum(x[, j]^2))
    if (scale > 0) x[, j] <- x[, j] / scale
  }
  w <- ifelse(boundary, 1, weight)
  root <- sqrt(w)
  absorb <- function(v) {
    fixest::demean(v, f = groups, weights = w, tol = 1e-12, iter = 10000L, as.matrix = TRUE, notes = FALSE)
  }
  basis <- if (ncol(x)) qr(absorb(x) * root, tol = 1e-10)
  residual <- function(target) {
    r <- drop(absorb(matrix(target)))
    if (is.null(basis)) r else qr.resid(basis, r * root) / root
  }
  function(u, tol) {
    target <- u
    for (iter in seq_len(maxiter)) {
      fit <- target - residual(target)
      gap <- u[!boundary] - fit[!boundary]
      if (max(abs(gap), 0) < tol) return(u - fit)
      target[!boundary] <- target[!boundary] + gap
    }
    stop("Constrained least squares in the separation check did not converge; no result returned.")
  }
}

relu_iterate <- function(boundary, project, tol = 1e-4, maxiter = 1000L) {
  u <- as.numeric(boundary)
  delta <- tol + max(sum(u^2) / ceiling(.Machine$double.eps^(-0.5)), 1e-8)
  for (iter in seq_len(maxiter)) {
    residual <- project(u, 0.1 * delta)
    z <- u - residual
    if (max(abs(z[!boundary])) > delta) {
      stop("ReLU equality constraints failed; separation was not certified.")
    }
    z[!boundary] <- 0
    z[z >= -0.1 * delta & z <= delta] <- 0
    if (all(z[boundary] >= 0)) return(boundary & z > 0)
    residual[abs(residual) < 1e-8] <- 0
    if (min(residual[boundary]) >= 0) {
      z[boundary & residual > delta] <- 0
      return(boundary & z > 0)
    }
    u[boundary] <- pmax(z[boundary], 0)
  }
  stop("ReLU separation check reached its iteration limit; no result returned.")
}
