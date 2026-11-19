# Regressors: treatment cells, alone and interacted with centered covariates,
# plus covariates by period and, if they vary within units, by cohort. One
# description serves the estimation formula and the sparse design matrices used
# for separation, variances and aggregation.

jw_covariates <- function(d, cvar) {
  if (is.null(cvar)) return(d)
  x <- stats::model.matrix(cvar, d)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  # Centering within cohort-period groups changes how cell coefficients are
  # reported, not the fit.
  centered <- fixest::demean(x,
    f = list(fixest::to_integer(d$.jw_g, d$.jw_t)), weights = d$.jw_w,
    tol = 1e-10, iter = 10000L, as.matrix = TRUE, notes = FALSE
  )
  for (j in seq_len(ncol(x))) {
    d[[paste0(".jw_x", j)]] <- x[, j]
    d[[paste0(".jw_xc", j)]] <- centered[, j]
  }
  d
}

jw_design <- function(d, gmax) {
  k <- length(grep("^\\.jw_x[0-9]+$", names(d)))
  varying <- vapply(seq_len(k), function(j) {
    spread <- tapply(d[[paste0(".jw_x", j)]], d$.jw_id, function(z) max(z) - min(z))
    any(spread > 1e-12)
  }, logical(1))
  cohorts <- sort(unique(d$.jw_g[d$.jw_g > 0 & d$.jw_g < gmax]))
  tmin <- min(d$.jw_t)
  gmin <- min(d$.jw_g)

  cells <- paste0(".jw_cell::", levels(d$.jw_cell)[-1L])
  terms <- "i(.jw_cell, ref = 'base')"
  columns <- cells
  for (j in seq_len(k)) {
    x <- paste0(".jw_x", j)
    terms <- c(terms, sprintf("i(.jw_cell, .jw_xc%d, ref = 'base')", j), x,
      sprintf("i(.jw_t, %s, ref = %s)", x, tmin),
      if (varying[j]) {
        sprintf("i(.jw_g, %s, ref = %s, keep = c(%s))", x, gmin, toString(cohorts))
      }
    )
    columns <- c(columns, paste0(cells, ":.jw_xc", j), x,
      paste0(".jw_t::", setdiff(sort(unique(d$.jw_t)), tmin), ":", x),
      if (varying[j]) paste0(".jw_g::", setdiff(cohorts, gmin), ":", x)
    )
  }
  fixed <- c(".jw_id", ".jw_t", grep("^\\.jw_fe[0-9]+$", names(d), value = TRUE))
  list(
    formula = stats::as.formula(
      paste(".jw_y ~", paste(terms, collapse = " + "), "|", paste(fixed, collapse = " + ")),
      env = baseenv()
    ),
    columns = columns, covariates = k, varying = varying
  )
}

# Sparse columns of the design matrix, aligned with `names`: the treatment
# cells and their covariate interactions, or the remaining covariate terms.
jw_columns <- function(d, names, design, part = c("treatment", "static")) {
  part <- match.arg(part)
  entries <- if (part == "treatment") {
    rows <- which(d$.jw_cell != "base")
    cell <- paste0(".jw_cell::", d$.jw_cell[rows])
    c(
      list(entry(rows, cell, 1)),
      lapply(seq_len(design$covariates), function(j) {
        entry(rows, paste0(cell, ":.jw_xc", j), d[[paste0(".jw_xc", j)]][rows])
      })
    )
  } else {
    rows <- seq_len(nrow(d))
    unlist(lapply(seq_len(design$covariates), function(j) {
      x <- paste0(".jw_x", j)
      list(
        entry(rows, x, d[[x]]),
        entry(rows, paste0(".jw_t::", d$.jw_t, ":", x), d[[x]]),
        if (design$varying[j]) entry(rows, paste0(".jw_g::", d$.jw_g, ":", x), d[[x]])
      )
    }), recursive = FALSE)
  }
  i <- as.integer(unlist(lapply(entries, `[[`, "i")))
  j <- match(unlist(lapply(entries, `[[`, "name")), names)
  x <- as.numeric(unlist(lapply(entries, `[[`, "x")))
  found <- !is.na(j)
  Matrix::sparseMatrix(
    i = i[found], j = j[found], x = x[found],
    dims = c(nrow(d), length(names)), dimnames = list(NULL, names)
  )
}

jw_design_matrix <- function(d, names, design) {
  jw_columns(d, names, design, "static") + jw_columns(d, names, design, "treatment")
}

entry <- function(i, name, x) list(i = i, name = rep_len(name, length(i)), x = rep_len(x, length(i)))
