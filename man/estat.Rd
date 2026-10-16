\name{estat}
\alias{estat}
\alias{print.jwdid_estat}
\title{Aggregate Average Treatment Effects}
\description{
Aggregates counterfactual prediction differences from a \code{jwdid} model
and computes clustered delta-method standard errors.
}
\usage{
estat(object,
  type = c("simple", "group", "calendar", "event", "attgt", "any"),
  predict = c("xb", "mu"), weights = NULL,
  subset = NULL, over = NULL, window = NULL, cwindow = NULL,
  asis = FALSE, level = 0.95, pretrend = FALSE)
}
\arguments{
\item{object}{A model fitted with \code{jwdid()}.}
\item{type}{Aggregate effects in periods after the treatment baseline overall (\code{simple}), by
cohort (\code{group}), or by period (\code{calendar}); aggregate by event time
(\code{event}) or cohort-period (\code{attgt}); or pool treated-cohort
observations before and after treatment (\code{any}). Event aggregation includes
pre-treatment effects when \code{never = TRUE}.}
\item{predict}{Prediction scale. The default, \code{xb}, returns log ratios
of conditional means. \code{mu} returns differences in outcome units.}
\item{weights}{Positive aggregation weights. Defaults to the estimation weights.
Accepts a vector aligned with the retained estimation data, a column name or
one-sided formula evaluated on that data.}
\item{subset}{Logical vector aligned with retained estimation data, column name,
or one-sided formula selecting observations for aggregation.}
\item{over}{Vector, column name or one-sided formula defining additional groups.}
\item{window}{Two endpoints restricting event time, in time-variable units.}
\item{cwindow}{Two endpoints collapsing event-time tails into the endpoints.
Use either \code{window} or \code{cwindow}.}
\item{asis}{Use the observed treatment indicator for the treated prediction
instead of setting it to one.}
\item{level}{Confidence level between zero and one.}
\item{pretrend}{For event aggregation with \code{never = TRUE}, test the joint
null that effects before the baseline are zero.}
}
\details{
Effects compare predictions with the treatment indicator set to zero and one,
holding covariates, offsets and fixed effects fixed. With \code{asis = TRUE},
the second prediction uses the observed treatment indicator.
The observation-level log ratios (\code{xb}) or outcome differences
(\code{mu}) are averaged using the aggregation weights.
Standard errors use the clustering specified when fitting the model.
Confidence intervals and p-values use the normal approximation.

With \code{never = TRUE}, event aggregation returns the baseline as a zero
effect with zero standard error. Its statistic and p-value are undefined.
A pretrend test with collapsed windows uses the reported collapsed event groups.
}
\value{
A data frame of class \code{jwdid_estat} containing group identifiers,
\code{estimate}, \code{std.error}, \code{statistic}, \code{p.value},
\code{conf.low}, \code{conf.high}, \code{n} and \code{weight_sum}.
The printed table shows group identifiers, estimates, standard errors and
p-values; all columns remain available in the returned data frame.
Attributes contain \code{vcov}, \code{gradient}, \code{type}, \code{predict}
and \code{df}. If requested, \code{pretrend} stores the joint chi-squared test.
}
\seealso{\code{\link{jwdid}}, \code{\link{plot.jwdid_estat}}}
\examples{
data(gravity)
fit <- jwdid(trade, ivar = exporter^importer, tvar = year, gvar = first_year,
  fevar = ~ exporter^year + importer^year, never = TRUE, data = gravity)
es <- estat(fit, "event", pretrend = TRUE)
print(attr(es, "pretrend"))
es_window <- estat(fit, "event", window = c(-2, 3))
plot(es_window)
}
