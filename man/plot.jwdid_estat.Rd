\name{plot.jwdid_estat}
\alias{plot.jwdid_estat}
\title{Plot Aggregated Treatment Effects}
\description{
Plots event, calendar or cohort effects and confidence intervals with
\code{fixest::coefplot()}.
}
\usage{
\method{plot}{jwdid_estat}(x, xlab = NULL,
  ylab = "Estimated treatment effect", main = "", ...)
}
\arguments{
\item{x}{A \code{jwdid_estat} table from event, calendar or group aggregation.}
\item{xlab}{Horizontal axis label. Defaults to the grouping column name.}
\item{ylab}{Vertical axis label.}
\item{main}{Plot title.}
\item{...}{Additional arguments to \code{fixest::coefplot()}.}
}
\details{Tables aggregated with \code{over} must first be reduced to one group,
with the \code{over} column removed.}
\value{Returns \code{x} invisibly after drawing the plot.}
\seealso{\code{\link{estat}}, \code{\link[fixest]{coefplot}}}
\examples{
data(gravity)
fit <- jwdid(trade, ivar = exporter^importer, tvar = year, gvar = first_year,
  fevar = ~ exporter^year + importer^year, never = TRUE, data = gravity)
es <- estat(fit, "event")
plot(es)
}
