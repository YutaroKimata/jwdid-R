.jwdid_vector <- function(x, data, default = NULL) {
  if (is.null(x)) {
    return(default)
  }
  if (inherits(x, "formula")) x <- eval(x[[2]], data, environment(x))
  if (is.character(x) && length(x) == 1L && x %in% names(data)) x <- data[[x]]
  if (length(x) != nrow(data)) stop("Vector length must match nrow(data).")

  x
}

.jwdid_design <- function(data, coefficients, treatment_only = FALSE) {
  blocks <- list()
  add <- function(columns, index, value, treatment) {
    column <- match(columns, coefficients, nomatch = 0L)[index]
    row <- which(column > 0L)
    if (length(row)) {
      blocks[[length(blocks) + 1L]] <<- list(
        row = row, column = column[row], value = value[row], treatment = treatment
      )
    }
  }

  cell <- as.integer(data$.jw_cell)
  cell_names <- paste0("fixest::.jw_cell::", levels(data$.jw_cell), ":.jw_tr")
  add(cell_names, cell, rep(1, nrow(data)), TRUE)
  for (name in grep("^\\.jw_x[[:digit:]]+$", names(data), value = TRUE)) {
    add(paste0(cell_names, ":", name), cell, data[[name]], TRUE)
  }

  if (!treatment_only) {
    times <- unique(data$.jw_t)
    cohorts <- unique(data$.jw_g)
    time <- match(data$.jw_t, times)
    cohort <- match(data$.jw_g, cohorts)
    for (name in grep("^\\.jw_raw[[:digit:]]+$", names(data), value = TRUE)) {
      add(name, rep.int(1L, nrow(data)), data[[name]], FALSE)
      add(paste0("fixest::.jw_t::", times, ":", name), time, data[[name]], FALSE)
      add(paste0("fixest::.jw_g::", cohorts, ":", name), cohort, data[[name]], FALSE)
    }
  }

  blocks
}

jwdid <- function(
  yvar, ivar, tvar, gvar, cvar = NULL, fevar = NULL,
  never = FALSE, anticipation = 1L,
  hettype = c("timecohort", "time", "cohort", "event", "twfe"),
  data, subset = NULL, weights = NULL, offset = NULL,
  cluster = NULL, ...
) {
  outcome_expr <- substitute(yvar)
  if (is.character(outcome_expr)) outcome_expr <- str2lang(outcome_expr)
  fml <- as.formula(
    call("~", outcome_expr, if (is.null(cvar)) 1 else cvar[[2]]),
    env = if (is.null(cvar)) parent.frame() else environment(cvar)
  )
  hettype <- match.arg(hettype)
  tvar <- as.character(substitute(tvar))
  gvar <- as.character(substitute(gvar))
  unit <- substitute(ivar)
  if (is.character(unit)) unit <- str2lang(unit)
  ivar <- if (is.symbol(unit)) as.character(unit) else ".jw_id"

  split_plus <- function(x) {
    if (is.call(x) && identical(x[[1]], as.name("+"))) {
      return(c(split_plus(x[[2]]), split_plus(x[[3]])))
    }
    list(x)
  }

  if (!is.numeric(anticipation) || length(anticipation) != 1L ||
    is.na(anticipation) || anticipation < 1 || anticipation != floor(anticipation)) {
    stop("anticipation is a positive integer; 1 means the usual g-1 baseline.")
  }

  d <- as.data.frame(data)
  if (any(startsWith(names(d), ".jw_"))) stop("Column names beginning .jw_ are reserved.")
  if (identical(ivar, ".jw_id")) {
    d$.jw_id <- fixest::to_integer(inputs = unname(as.list(d[all.vars(unit)])))
  }

  w <- .jwdid_vector(weights, d, rep(1, nrow(d)))
  off <- .jwdid_vector(offset, d, rep(0, nrow(d)))
  keep <- .jwdid_vector(subset, d, rep(TRUE, nrow(d)))
  if (any(w <= 0, na.rm = TRUE)) stop("Weights must be strictly positive.")
  vars <- unique(c(all.vars(fml), ivar, tvar, gvar, all.vars(fevar), all.vars(cluster)))
  keep <- !is.na(keep) & keep & stats::complete.cases(d[vars], w) & is.finite(w)
  source_rows <- which(keep)
  if (!all(keep)) {
    d <- d[keep, , drop = FALSE]
    w <- w[keep]
    off <- off[keep]
  }
  if (!is.numeric(d[[tvar]]) || any(!is.finite(d[[tvar]]))) {
    stop("The panel time variable must contain finite numeric periods.")
  }
  if (anyDuplicated(fixest::to_integer(d[[ivar]], d[[tvar]]))) {
    stop("Duplicate individual-time observations in panel.")
  }

  d$.jw_g <- d[[gvar]]
  if (!is.numeric(d$.jw_g) || any(!is.finite(d$.jw_g)) || any(d$.jw_g < 0)) {
    stop("gvar must be numeric, nonnegative; 0 denotes never treated.")
  }
  if (any(d$.jw_g != d$.jw_g[match(d[[ivar]], d[[ivar]])])) {
    stop("gvar must be constant within individuals (absorbing treatment).")
  }

  d$.jw_t <- d[[tvar]]
  time_order <- order(d$.jw_t)
  min_t <- d$.jw_t[time_order[match(d[[ivar]], d[[ivar]][time_order])]]
  keep <- !(d$.jw_g > 0 & min_t >= d$.jw_g)
  gmax <- if (any(d$.jw_g[keep] == 0)) Inf else max(d$.jw_g[keep])
  keep <- keep & d$.jw_t < gmax
  periods <- sort(unique(c(d$.jw_g[keep], d$.jw_t[keep])))
  periods <- periods[periods > 0]
  if (length(periods) < 2L) stop("At least two usable positive periods are required.")

  gap <- min(diff(periods))
  base <- d$.jw_g - gap * anticipation
  keep <- keep & !(d$.jw_g > 0 & min_t >= d$.jw_g - gap * (anticipation - 1L))
  if (!all(keep)) {
    source_rows <- source_rows[keep]
    d <- d[keep, , drop = FALSE]
    w <- w[keep]
    base <- base[keep]
    off <- off[keep]
  }
  d$.jw_tr <- as.numeric(
    d$.jw_g < gmax & if (never) TRUE else d$.jw_g > 0 & d$.jw_t >= base
  )
  d$.jw_etr <- d$.jw_g > 0 & d$.jw_t > base
  d$.jw_event <- d$.jw_t - d$.jw_g
  eligible <- d$.jw_g > 0 & d$.jw_g < gmax & if (never) d$.jw_t != base else d$.jw_t > base
  post <- if (hettype %in% c("time", "cohort", "twfe")) {
    ifelse(d$.jw_t < base, "pre", "post")
  }
  cell <- switch(hettype,
    timecohort = paste0("g", d$.jw_g, "_t", d$.jw_t),
    time = paste0(post, "_t", d$.jw_t),
    cohort = paste0("g", d$.jw_g, "_", post),
    event = paste0("e", d$.jw_event),
    twfe = post
  )
  cell[!eligible] <- "base"
  d$.jw_cell <- factor(cell, levels = c("base", sort(unique(cell[cell != "base"]))))
  if (nlevels(d$.jw_cell) < 2L) stop("No identified treatment cells.")

  treatment <- "fixest::i(.jw_cell, .jw_tr, ref = 'base')"
  rhs <- treatment
  if (!is.null(cvar)) {
    x <- stats::model.matrix(stats::delete.response(stats::terms(fml)), data = d)
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
    if (ncol(x)) {
      center_groups <- switch(hettype,
        timecohort = fixest::to_integer(d$.jw_g, d$.jw_t),
        time = fixest::to_integer(d$.jw_t, ifelse(eligible, post, "base")),
        cohort = fixest::to_integer(d$.jw_g, ifelse(eligible, post, "base")),
        event = d$.jw_event,
        twfe = ifelse(eligible, post, "base")
      )
      centered <- fixest::demean(
        x,
        f = list(center_groups), weights = w,
        tol = 1e-10, iter = 10000L, as.matrix = TRUE, notes = FALSE
      )
      units <- as.factor(d[[ivar]])
      for (j in seq_len(ncol(x))) {
        nm <- sprintf(".jw_x%03d", j)
        d[[nm]] <- centered[, j]
        rhs <- c(rhs, paste0(treatment, ":", nm))
        raw <- sprintf(".jw_raw%03d", j)
        d[[raw]] <- x[, j]
        rhs <- c(rhs, raw, paste0("fixest::i(.jw_t, ", raw, ", ref = ", min(d$.jw_t), ")"))
        if (any(vapply(
          split(x[, j], units),
          function(z) max(z) - min(z), numeric(1)
        ) > 1e-12)) {
          rhs <- c(rhs, paste0(
            "fixest::i(.jw_g, ", raw,
            ", ref = ", min(d$.jw_g), ", keep = c(",
            paste(sort(unique(d$.jw_g[d$.jw_g > 0 & d$.jw_g < gmax])), collapse = ","), "))"
          ))
        }
      }
      rm(centered, units, center_groups)
    }
    rm(x)
  }

  if (is.null(cluster)) cluster <- stats::reformulate(ivar, env = environment(fml))
  linear <- as.formula(
    paste(deparse1(fml[[2]], collapse = ""), "~", paste(rhs, collapse = " + ")),
    env = environment(fml)
  )
  absorb <- c(ivar, ".jw_t", if (!is.null(fevar)) deparse1(fevar[[2]], collapse = ""))
  full <- stats::as.formula(
    paste(deparse1(linear, collapse = ""), "|", paste(absorb, collapse = " + ")),
    env = environment(linear)
  )
  mf <- stats::model.frame(fml, d, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (!is.numeric(y) || any(y < 0, na.rm = TRUE)) {
    stop("PPML requires nonnegative numeric outcomes.")
  }
  required <- unique(c(all.vars(full), all.vars(cluster)))
  keep <- stats::complete.cases(d[required], mf, w, off) & is.finite(y) & is.finite(off)
  if (!all(keep)) {
    source_rows <- source_rows[keep]
    d <- d[keep, , drop = FALSE]
    w <- w[keep]
    off <- off[keep]
    y <- y[keep]
  }
  rm(mf, time_order, min_t, base, post, cell, eligible, keep)

  separated_rows <- integer()
  group_vectors <- function(fml, data) {
    expr <- if (inherits(fml, "formula")) {
      fml[[length(fml)]]
    } else {
      str2lang(paste(fml, collapse = " + "))
    }
    lapply(split_plus(expr), function(term) {
      if (grepl("[\\[\\]]", deparse1(term, collapse = ""))) {
        stop("Varying-slope fixed effects are not supported in this release.")
      }
      if (is.call(term) && identical(term[[1]], as.name("^"))) {
        fixest::to_integer(inputs = unname(as.list(data[all.vars(term)])))
      } else {
        fixest::to_integer(eval(term, data, parent.frame()))
      }
    })
  }
  groups <- group_vectors(absorb, d)
  relu_separation <- function(y, data, groups, tol = 1e-4, maxiter = 100L) {
    boundary <- y == 0
    w <- ifelse(boundary, 1, ceiling(.Machine$double.eps^(-0.5)))
    u <- as.numeric(boundary)
    qr_x <- NULL
    if (is.null(cvar)) {
      groups <- c(groups, list(fixest::to_integer(data$.jw_cell)))
    } else {
      sqrt_w <- sqrt(w)
      x <- stats::model.matrix(linear, data = data)
      x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
      for (j in seq_len(ncol(x))) {
        scale <- sqrt(sum(x[, j]^2))
        if (scale > 0) x[, j] <- x[, j] / scale
      }
      x <- if (ncol(x)) {
        fixest::demean(
          x,
          f = groups, weights = w,
          tol = 1e-12, iter = 10000L, as.matrix = TRUE, notes = FALSE
        )
      } else {
        x
      }
      if (ncol(x)) {
        x <- x * sqrt_w
        qr_x <- qr(x, tol = 1e-10)
      }
      rm(x)
    }
    epsilon <- max(sum(u^2) / max(w), 1e-8)
    delta <- tol + epsilon
    for (iter in seq_len(maxiter)) {
      ud <- drop(fixest::demean(
        matrix(u),
        f = groups, weights = w,
        tol = 1e-12, iter = 10000L, as.matrix = TRUE, notes = FALSE
      ))
      residual <- if (is.null(qr_x)) ud else qr.resid(qr_x, ud * sqrt_w) / sqrt_w
      z <- u - residual
      if (max(abs(z[!boundary])) > delta) {
        stop("ReLU equality constraints failed; separation was not certified.")
      }
      z[!boundary] <- 0
      z[z >= -0.1 * delta & z <= delta] <- 0
      if (all(z[boundary] >= 0)) {
        return(boundary & z > 0)
      }
      residual[abs(residual) < 1e-8] <- 0
      if (min(residual[boundary]) >= 0) {
        z[boundary & residual > delta] <- 0
        return(boundary & z > 0)
      }
      u[boundary] <- pmax(z[boundary], 0)
    }
    stop("ReLU separation check reached its iteration limit; no result returned.")
  }

  repeat {
    bad <- rep(FALSE, nrow(d))
    for (g in groups) bad <- bad | rowsum(y, g)[g, 1L] == 0
    if (!any(bad) && any(y == 0)) bad <- relu_separation(y, d, groups)
    if (!any(bad)) break
    separated_rows <- c(separated_rows, source_rows[bad])
    if (all(bad)) stop("No observations remain after separation removal.")
    source_rows <- source_rows[!bad]
    d <- d[!bad, , drop = FALSE]
    w <- w[!bad]
    off <- off[!bad]
    y <- y[!bad]
    groups <- lapply(groups, function(g) fixest::to_integer(g[!bad]))
  }
  rm(groups, y)

  cg <- group_vectors(cluster, d)
  fit <- do.call(
    fixest::fepois,
    c(
      list(
        fml = full, data = d,
        weights = if (is.null(weights)) NULL else w,
        offset = if (is.null(offset)) NULL else off,
        vcov = if (length(cg) == 1L) {
          cluster
        } else {
          fixest::vcov_cluster(cluster = cluster, vcov_fix = FALSE)
        },
        ssc = fixest::ssc(K.adj = FALSE, G.adj = TRUE, G.df = "min"),
        fixef.rm = "none", data.save = TRUE
      ),
      utils::modifyList(
        list(glm.tol = 1e-10, glm.iter = 10000L, fixef.tol = 1e-9),
        list(...)
      )
    )
  )
  parent.env(fit$call_env) <- environment(fml)
  if (!isTRUE(fit$convStatus)) stop("PPML did not converge; no result returned.")

  used <- fixest::obs(fit)
  if (length(used) != nrow(d)) {
    d <- d[used, , drop = FALSE]
    w <- w[used]
    source_rows <- source_rows[used]
    cg <- lapply(cg, function(g) fixest::to_integer(g[used]))
  }

  b <- stats::coef(fit)
  v <- fit$cov.scaled[, , drop = FALSE]
  q <- stats::fitted(fit) * w
  q_sum <- sum(q)
  cell_coefficients <- NULL
  if (is.null(cvar)) {
    cell_coefficients <- match(
      paste0("fixest::.jw_cell::", levels(d$.jw_cell), ":.jw_tr"), names(b)
    )
    cell <- cell_coefficients[as.integer(d$.jw_cell)]
    active <- !is.na(cell)
    sums <- rowsum(q[active] * d$.jw_tr[active], cell[active])
    xm <- numeric(length(b))
    xm[as.integer(rownames(sums))] <- sums[, 1L] / q_sum
    if (length(cg) > 1L) {
      counts <- tabulate(cell, nbins = length(b))
      scales <- sqrt(counts * (1 - counts / nrow(d)) / (nrow(d) - 1L))
    }
  } else {
    xm <- numeric(length(b))
    if (length(cg) > 1L) scales <- numeric(length(b))
    design <- .jwdid_design(d, names(b))
    for (block in design) {
      value <- block$value
      if (block$treatment) value <- value * d$.jw_tr[block$row]
      sums <- rowsum(value * q[block$row], block$column)
      columns <- as.integer(rownames(sums))
      xm[columns] <- sums[, 1L] / q_sum
      if (length(cg) > 1L) {
        means <- numeric(length(b))
        means[columns] <- rowsum(value, block$column)[, 1L] / nrow(d)
        squares <- rowsum((value - means[block$column])^2, block$column)[, 1L]
        counts <- tabulate(block$column, length(b))[columns]
        scales[columns] <- sqrt(
          (squares + (nrow(d) - counts) * means[columns]^2) / (nrow(d) - 1L)
        )
      }
    }
    rm(design, block, value)
  }
  rm(q)

  residual <- stats::residuals(fit, type = "response") * w
  scores <- fixest::estfun(fit)
  bread <- fixest::bread(fit) / stats::nobs(fit)
  sizes <- vapply(cg, max, integer(1))
  score_cross <- numeric(length(b))
  residual_variance <- 0
  for (k in seq_along(cg)) {
    for (idx in utils::combn(seq_along(cg), k, simplify = FALSE)) {
      g <- if (k == 1L) cg[[idx]] else fixest::to_integer(inputs = cg[idx])
      score_sum <- rowsum(scores, g, reorder = FALSE)
      residual_sum <- drop(rowsum(residual, g, reorder = FALSE)) / q_sum
      score_cross <- score_cross + (-1)^(k + 1L) * drop(crossprod(score_sum, residual_sum))
      residual_variance <- residual_variance + (-1)^(k + 1L) * sum(residual_sum^2)
    }
  }
  adjustment <- min(sizes) / (min(sizes) - 1L)
  cross <- drop(crossprod(bread, score_cross)) * adjustment
  vx <- drop(v %*% xm)
  intercept_vcov <- c(
    cross - vx,
    `_cons` = residual_variance * adjustment - 2 * sum(xm * cross) + sum(xm * vx)
  )
  if (length(cg) > 1L) {
    v <- rbind(cbind(v, `_cons` = intercept_vcov[seq_along(b)]), `_cons` = intercept_vcov)
    scales <- c(pmax(scales, 1e-3), 1)
    scaled_v <- v * outer(scales, scales)
    eig <- eigen(scaled_v, symmetric = TRUE)
    if (min(eig$values) < -1e-10 * max(abs(eig$values))) {
      warning("Multiway covariance was not positive semidefinite; eigenvalues clipped.")
      v <- (eig$vectors %*% (pmax(eig$values, 0) * t(eig$vectors))) / outer(scales, scales)
    }
    dimnames(v) <- list(c(names(b), "_cons"), c(names(b), "_cons"))
    intercept_vcov <- v[nrow(v), ]
    v <- v[-nrow(v), -ncol(v), drop = FALSE]
  }

  fit$cov.scaled <- v
  fit$se <- sqrt(diag(v))
  fit$coeftable[, 2L] <- fit$se
  fit$coeftable[, 3L] <- b / fit$se
  fit$coeftable[, 4L] <- 2 * stats::pnorm(-abs(fit$coeftable[, 3L]))
  fit$summary_flags$vcov <- cluster

  fit$jwdid <- list(
    data = d, weights = w,
    source_rows = source_rows, ivar = unit, tvar = tvar, gvar = gvar,
    cvar = cvar, fevar = fevar, never = never, gap = gap, anticipation = anticipation,
    hettype = hettype, cluster = cluster, generated_formula = stats::formula(fit),
    cell_coefficients = cell_coefficients, intercept_vcov = intercept_vcov,
    separated_rows = sort(unique(separated_rows))
  )

  fit
}

estat <- function(
  object, type = c("simple", "group", "calendar", "event", "attgt", "any"),
  predict = c("xb", "mu"), weights = NULL,
  subset = NULL, over = NULL, window = NULL, cwindow = NULL,
  asis = FALSE, level = 0.95, pretrend = FALSE
) {
  type <- match.arg(type)
  predict <- match.arg(predict)
  if (length(level) != 1L || level <= 0 || level >= 1) stop("level must be in (0,1).")

  meta <- object$jwdid
  d <- meta$data
  w <- .jwdid_vector(weights, d, meta$weights)
  if (any(!is.finite(w) | w <= 0)) stop("Aggregation weights must be finite and positive.")
  keep <- .jwdid_vector(subset, d, rep(TRUE, nrow(d)))
  keep <- !is.na(keep) & keep
  if (type %in% c("simple", "group", "calendar")) keep <- keep & d$.jw_etr
  if (type == "event") {
    keep <- keep & if (meta$never) d$.jw_g > 0 & d$.jw_tr != 0 else d$.jw_etr
  }
  if (type %in% c("attgt", "any")) keep <- keep & d$.jw_g > 0
  ev <- d$.jw_event
  if (!is.null(window) && !is.null(cwindow)) stop("Use either window or cwindow.")
  bounds <- if (!is.null(window)) window else cwindow
  if (!is.null(bounds) && (length(bounds) != 2L || bounds[1L] > bounds[2L])) {
    stop("A window is c(lower, upper) in calendar-time units.")
  }
  if (!is.null(window)) keep <- keep & ev >= window[1L] & ev <= window[2L]
  if (!is.null(cwindow)) ev <- pmin(pmax(ev, cwindow[1L]), cwindow[2L])

  keys <- switch(type,
    simple = data.frame(aggregate = rep("simple", nrow(d))),
    any = data.frame(aggregate = rep("any", nrow(d))),
    group = data.frame(cohort = d$.jw_g),
    calendar = data.frame(time = d$.jw_t),
    event = data.frame(event = ev),
    attgt = data.frame(cohort = d$.jw_g, time = d$.jw_t)
  )
  if (!is.null(over)) {
    ov <- .jwdid_vector(over, d)
    keep <- keep & !is.na(ov)
    keys$over <- ov
  }
  if (!any(keep)) stop("No observations in the requested aggregation.")

  b <- stats::coef(object)
  key <- do.call(interaction, c(keys, list(drop = TRUE, lex.order = TRUE)))
  key <- droplevels(key[keep])
  group <- as.integer(key)
  first <- which(keep)[match(seq_len(nlevels(key)), group)]
  out <- keys[first, , drop = FALSE]
  n <- stats::setNames(tabulate(group, nlevels(key)), levels(key))
  w <- w[keep]
  weight_sum <- stats::setNames(drop(rowsum(w, group)), levels(key))
  vars <- all.vars(stats::formula(object, type = "linear")[[3L]])
  d <- d[keep, vars, drop = FALSE]
  mu <- if (predict == "mu") stats::fitted(object)[keep] else NULL

  if (predict == "xb" || is.null(meta$cvar)) {
    cell <- fixest::to_integer(group, d$.jw_cell)
    first <- match(seq_len(max(cell)), cell)
    cell_weight <- drop(rowsum(w, cell))
    compressed <- d[first, , drop = FALSE]
    for (name in grep("^\\.jw_x[[:digit:]]+$", names(d), value = TRUE)) {
      compressed[[name]] <- drop(rowsum(d[[name]] * w, cell)) / cell_weight
    }
    if (predict == "mu") mu <- drop(rowsum(mu * w, cell)) / cell_weight
    d <- compressed
    w <- cell_weight
    group <- group[first]
  }

  g <- matrix(
    0, nlevels(key), length(b) + 1L,
    dimnames = list(levels(key), c(names(b), "_cons"))
  )
  if (is.null(meta$cvar)) {
    columns <- meta$cell_coefficients[as.integer(d$.jw_cell)]
    active <- !is.na(columns)
    target <- if (asis) d$.jw_tr else 1
    if (predict == "xb") {
      weighted_gradient <- w * target
    } else {
      delta <- numeric(nrow(d))
      delta[active] <- b[columns[active]]
      eta <- log(mu)
      p0 <- exp(eta - d$.jw_tr * delta)
      p1 <- exp(eta + (target - d$.jw_tr) * delta)
      weighted_effects <- (p1 - p0) * w
      weighted_gradient <- target * p1 * w
      g[, length(b) + 1L] <- drop(rowsum(weighted_effects, group))
    }
    g[cbind(group[active], columns[active])] <- weighted_gradient[active]
  } else {
    design <- .jwdid_design(d, names(b), treatment_only = predict == "xb")
    target <- if (asis) d$.jw_tr else 1
    if (predict == "xb") {
      treatment_weight <- w * target
    } else {
      delta <- numeric(nrow(d))
      for (block in design) {
        if (block$treatment) {
          delta[block$row] <- delta[block$row] + block$value * b[block$column]
        }
      }
      eta <- log(mu)
      p0 <- exp(eta - d$.jw_tr * delta)
      p1 <- exp(eta + (target - d$.jw_tr) * delta)
      effect_weight <- (p1 - p0) * w
      treatment_weight <- target * p1 * w
      g[, length(b) + 1L] <- drop(rowsum(effect_weight, group))
    }
    for (block in design) {
      weight <- if (block$treatment) treatment_weight else effect_weight
      index <- group[block$row] + (block$column - 1L) * nrow(g)
      sums <- rowsum(block$value * weight[block$row], index, reorder = FALSE)
      index <- as.integer(rownames(sums))
      g[index] <- g[index] + sums[, 1L]
    }
  }
  g <- g / weight_sum
  slopes <- g[, seq_along(b), drop = FALSE]
  estimates <- if (predict == "xb") drop(slopes %*% b) else g[, ncol(g)]
  av <- tcrossprod(slopes %*% stats::vcov(object), slopes)
  if (predict == "mu") {
    intercept <- g[, length(b) + 1L]
    cross <- drop(slopes %*% meta$intercept_vcov[seq_along(b)])
    av <- av + tcrossprod(cross, intercept) + tcrossprod(intercept, cross) +
      meta$intercept_vcov[length(b) + 1L] * tcrossprod(intercept)
  }
  se <- sqrt(pmax(diag(av), 0))
  critical <- stats::qnorm((1 + level) / 2)
  statistic <- estimates / se
  pvalue <- 2 * stats::pnorm(-abs(statistic))
  rownames(out) <- NULL
  out[c(
    "estimate", "std.error", "statistic", "p.value",
    "conf.low", "conf.high", "n", "weight_sum"
  )] <- list(
    estimates, se, statistic, pvalue, estimates - critical * se, estimates + critical * se,
    n, weight_sum
  )
  out <- structure(out, vcov = av, gradient = g, type = type, predict = predict, df = Inf)
  if (pretrend) {
    if (type != "event" || !meta$never) {
      stop("pretrend requires event aggregation and never=TRUE.")
    }
    idx <- which(out$event < -meta$gap * meta$anticipation)
    if (!length(idx)) stop("No pre-treatment effects in the requested window.")
    vv <- av[idx, idx, drop = FALSE]
    ee <- estimates[idx]
    eig <- eigen(vv, symmetric = TRUE)
    retain <- eig$values > max(eig$values, 0) * 1e-9
    rank <- sum(retain)
    if (!rank) stop("No estimable pre-treatment effects for the joint test.")
    chi2 <- sum(
      drop(crossprod(eig$vectors[, retain, drop = FALSE], ee))^2 / eig$values[retain]
    )
    attr(out, "pretrend") <- list(
      statistic = chi2, df = rank, df_residual = Inf,
      p.value = stats::pchisq(chi2, rank, lower.tail = FALSE)
    )
  }

  class(out) <- c("jwdid_estat", "data.frame")
  out
}

print.jwdid_estat <- function(x, ...) {
  cols <- setdiff(names(x), c("statistic", "conf.low", "conf.high", "n", "weight_sum"))
  print.data.frame(x[cols], ...)
  invisible(x)
}

plot.jwdid_estat <- function(
  x, xlab = NULL, ylab = "Estimated treatment effect", main = "", ...
) {
  key <- intersect(c("event", "time", "cohort"), names(x))[1L]
  if ("over" %in% names(x)) stop("Subset to one over group before plotting.")
  ord <- order(x[[key]])
  fixest::coefplot(
    x$estimate[ord],
    x = x[[key]][ord],
    ci_low = x$conf.low[ord], ci_high = x$conf.high[ord],
    xlab = if (is.null(xlab)) key else xlab, ylab = ylab, main = main, ...
  )
  invisible(x)
}
