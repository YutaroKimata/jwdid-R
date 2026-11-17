# jwdid-R

An unofficial R implementation of Stata's
[jwdid](https://github.com/friosavila/stpackages/tree/738defceb5413face7051754b259c82798c67e52/jwdid),
estimated with [fixest](https://lrberge.github.io/fixest/).

## Install

```r
install.packages("remotes")
remotes::install_github("YutaroKimata/jwdid-R")
```

## Gravity example

```r
library(jwdid.r)

fit <- jwdid(
  yvar = trade,
  ivar = exporter^importer,
  tvar = year,
  gvar = first_year,
  fevar = ~ exporter^year + importer^year,
  never = TRUE,
  data = gravity
)
fit
```

This estimates the effects of regional trade agreements on bilateral trade in
the simulated `gravity` panel, using dyads that never adopt an agreement as
controls, with dyad, exporter-year and importer-year fixed effects.

## Aggregated effects

```r
estat(fit, "simple")
es <- estat(fit, "event")
plot(es)
```

`estat()` averages effects over treated observations overall, by cohort,
period, event time or cohort-period, on the log scale (`predict = "xb"`) or in
outcome units (`predict = "mu"`), with clustered delta-method standard errors.

See `?jwdid` and `?estat` for additional options.

[MIT license](LICENSE.md).
