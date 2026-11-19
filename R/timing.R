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

# Post-treatment indicator and cohort-period cells.
jw_cells <- function(t, g, gap, gmax, anticipation, never) {
  lag <- gap * anticipation
  lead <- gap * (anticipation - 1L)
  member <- g > 0 & g < gmax & if (never) t != g - lag else t >= g - lead
  label <- ifelse(member, paste0("g", g, "_t", t), "base")
  # Cells run by cohort and period, so the estimator drops collinear cells from the end.
  cells <- setdiff(unique(label[order(g, t)]), "base")
  list(
    .jw_etr = g > 0 & t > g - lag,
    .jw_cell = factor(label, c("base", cells))
  )
}
