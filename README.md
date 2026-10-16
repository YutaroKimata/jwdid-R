# jwdid-R

An unofficial R implementation of Stata's
[jwdid](https://github.com/friosavila/stpackages/tree/738defceb5413face7051754b259c82798c67e52/jwdid).

## Install

```r
install.packages("remotes")
remotes::install_github("YutaroKimata/jwdid-R")
```

## Gravity example

Using the included simulated trade panel:

```r
library(jwdid.r)

data(gravity)

fit <- jwdid(
  yvar = trade,
  ivar = exporter^importer,
  tvar = year,
  gvar = first_year,
  fevar = ~ exporter^year + importer^year,
  never = TRUE,
  data = gravity
)
```

This estimates the effects of regional trade agreements on bilateral trade,
using dyads that never adopt an agreement as controls. The model includes dyad,
exporter-year and importer-year fixed effects.

## ATT and event studies

```r
att <- estat(fit, "simple")
print(att)

es <- estat(fit, "event")
plot(es)
```

`att` summarizes the average effect on treated dyads. The event-study plot shows
effects before and after agreement adoption, with confidence intervals.

See `?jwdid` and `?estat` for additional options.

[MIT license](LICENSE.md).
