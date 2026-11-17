# Input handling: one data frame holds the user's columns plus internal `.jw_`
# columns, so every later subset keeps outcome, weights and groups aligned.

jw_panel <- function(data, roles, cvar, fevar, cluster, subset, weights, offset, env) {
  d <- as.data.frame(data)
  if (any(startsWith(names(d), ".jw_"))) stop("Column names beginning .jw_ are reserved.")
  d$.jw_row <- seq_len(nrow(d))
  d$.jw_y <- eval(roles$y, d, env)
  d$.jw_id <- group_codes(roles$id, d, env)
  d$.jw_t <- eval(roles$t, d, env)
  d$.jw_g <- eval(roles$g, d, env)
  d$.jw_w <- column(weights, d, rep(1, nrow(d)))
  d$.jw_off <- column(offset, d, rep(0, nrow(d)))
  fixed <- if (is.null(fevar)) list() else split_terms(fevar[[2]])
  clusters <- if (is.null(cluster)) list(quote(.jw_id)) else split_terms(cluster[[2]])
  for (k in seq_along(fixed)) d[[paste0(".jw_fe", k)]] <- group_codes(fixed[[k]], d, env)
  for (k in seq_along(clusters)) d[[paste0(".jw_cl", k)]] <- group_codes(clusters[[k]], d, env)

  if (any(d$.jw_w < 0, na.rm = TRUE)) stop("Weights must be nonnegative.")
  inputs <- c(".jw_y", ".jw_id", ".jw_t", ".jw_g", if (!is.null(cvar)) all.vars(cvar))
  keep <- column(subset, d, rep(TRUE, nrow(d)))
  keep <- !is.na(keep) & keep & stats::complete.cases(d[inputs]) & is.finite(d$.jw_w) & d$.jw_w > 0
  d <- d[keep, , drop = FALSE]

  if (!is.numeric(d$.jw_y)) stop("PPML requires a numeric outcome.")
  if (!is.numeric(d$.jw_t) || any(!is.finite(d$.jw_t))) {
    stop("The time variable must contain finite numeric periods.")
  }
  if (!is.numeric(d$.jw_g) || any(!is.finite(d$.jw_g) | d$.jw_g < 0)) {
    stop("gvar must be numeric and nonnegative; 0 marks never-treated units.")
  }
  if (anyDuplicated(d[c(".jw_id", ".jw_t")])) stop("Duplicate unit-time observations.")
  if (any(d$.jw_g != stats::ave(d$.jw_g, d$.jw_id, FUN = min))) {
    stop("gvar must be constant within units.")
  }
  d
}

# Integer codes for a grouping term; `a^b` combines variables. Missing inputs stay missing.
group_codes <- function(term, data, env) {
  values <- if (is.call(term) && identical(term[[1]], as.name("^"))) {
    unname(as.list(data[all.vars(term)]))
  } else {
    list(eval(term, data, env))
  }
  codes <- do.call(fixest::to_integer, values)
  codes[Reduce(`|`, lapply(values, is.na))] <- NA
  codes
}

split_terms <- function(expr) {
  if (grepl("[][]", deparse1(expr))) stop("Varying-slope fixed effects are not supported.")
  if (is.call(expr) && identical(expr[[1]], as.name("+"))) {
    return(c(split_terms(expr[[2]]), split_terms(expr[[3]])))
  }
  list(expr)
}

# A vector aligned with `data`, given as a vector, a column name or a one-sided formula.
column <- function(x, data, default) {
  if (is.null(x)) return(default)
  if (inherits(x, "formula")) x <- eval(x[[2]], data, environment(x))
  if (is.character(x) && length(x) == 1L && x %in% names(data)) x <- data[[x]]
  if (length(x) != nrow(data)) stop("Vectors must have one value per row of data.")
  x
}

as_expression <- function(x) if (is.character(x)) str2lang(x) else x
