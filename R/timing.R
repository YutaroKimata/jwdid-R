# Sample rules, as a mask over unit-period rows.
jw_timing <- function(id, t, g, anticipation) {
  first <- stats::ave(t, id, FUN = min)
  keep <- !(g > 0 & first >= g)
  gmax <- if (any(g[keep] == 0)) Inf else max(g[keep])
  keep <- keep & t < gmax
  periods <- sort(unique(c(g[keep], t[keep])))
  periods <- periods[periods > 0]
  if (length(periods) < 2L) stop("At least two usable positive periods are required.")
  gap <- min(diff(periods))
  keep <- keep & !(g > 0 & first >= g - gap * (anticipation - 1L))
  list(keep = keep, gap = gap, gmax = gmax)
}

# Treatment indicators and cells. `phase` is 1 before the baseline, 2 from the
# first anticipated period on and 0 otherwise.
jw_cells <- function(t, g, gap, gmax, anticipation, never, hettype) {
  lag <- gap * anticipation
  lead <- gap * (anticipation - 1L)
  phase <- ifelse(g > 0 & t < g - lag, 1L, ifelse(g > 0 & t >= g - lead, 2L, 0L))
  member <- g > 0 & switch(hettype,
    timecohort = g < gmax & if (never) t != g - lag else t >= g - lead,
    cohort = g < gmax & (phase == 2L | never & phase == 1L),
    event = if (never) t != g - lag else t > g - lag,
    phase == 2L | never & phase == 1L
  )
  side <- ifelse(phase == 1L, "pre", "post")
  label <- switch(hettype,
    timecohort = paste0("g", g, "_t", t),
    time = paste0(side, "_t", t),
    cohort = paste0("g", g, "_", side),
    event = paste0("e", t - g),
    twfe = side
  )
  label[!member] <- "base"
  # Cells run by cohort and period, so the estimator drops collinear cells from the end.
  sequence <- switch(hettype,
    timecohort = order(g, t), time = order(phase, t), cohort = order(g, phase),
    event = order(t - g), twfe = order(phase)
  )
  list(
    .jw_tr = as.numeric(g < gmax & (never | g > 0 & t >= g - lag)),
    .jw_etr = g > 0 & t > g - lag,
    .jw_phase = phase,
    .jw_cell = factor(label, c("base", setdiff(unique(label[sequence]), "base")))
  )
}
