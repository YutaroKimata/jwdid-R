set.seed(20261015)
countries <- 1:16
years <- 2000:2011
pairs <- t(combn(countries, 2))
first_year <- sample(rep(c(2003L, 2006L, 2009L, 0L), each = nrow(pairs) / 4))
dyads <- data.frame(
  exporter = c(pairs[, 1], pairs[, 2]),
  importer = c(pairs[, 2], pairs[, 1]),
  first_year = rep(first_year, 2),
  pair_fe = rep(rnorm(nrow(pairs), sd = 0.5), 2)
)
df <- merge(dyads, data.frame(year = years), by = NULL)
exporter_year <- matrix(rnorm(length(countries) * length(years), sd = 0.3), nrow = length(countries))
importer_year <- matrix(rnorm(length(countries) * length(years), sd = 0.3), nrow = length(countries))
period <- match(df$year, years)
b <- c(`2003` = 0.06, `2006` = 0.04, `2009` = 0.02, `0` = 0)
df$rta <- as.integer(df$first_year > 0 & df$year >= df$first_year)
df$tau <- unname(b[as.character(df$first_year)]) * pmax(df$year - df$first_year + 1, 0) * df$rta
df$mu <- exp(2 + df$pair_fe + 0.04 * (df$year - min(years)) +
  exporter_year[cbind(df$exporter, period)] + importer_year[cbind(df$importer, period)] + df$tau)
df$trade <- rpois(nrow(df), df$mu)
gravity <- df[order(df$exporter, df$importer, df$year),
  c("exporter", "importer", "year", "first_year", "rta", "trade", "tau", "mu")]
rownames(gravity) <- NULL
save(gravity, file = "data/gravity.rda", version = 2)
