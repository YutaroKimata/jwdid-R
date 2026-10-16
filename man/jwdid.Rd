\name{jwdid}
\alias{jwdid}
\title{Staggered Difference-in-Differences with PPML}
\description{
Estimates extended two-way fixed effects difference-in-differences models
using \code{fixest::fepois()}, with treatment cells constructed from a first
treatment period.
}
\usage{
jwdid(yvar, ivar, tvar, gvar, cvar = NULL, fevar = NULL,
  never = FALSE, anticipation = 1L,
  hettype = c("timecohort", "time", "cohort", "event", "twfe"),
  data, subset = NULL, weights = NULL, offset = NULL,
  cluster = NULL, ...)
}
\arguments{
\item{yvar}{Nonnegative numeric outcome, supplied as a column name or string.}
\item{ivar}{Unit identifier, supplied as a column name or string.
Combined identifiers can use \code{exporter^importer} or
\code{"exporter^importer"}.}
\item{tvar}{Numeric time variable, supplied as a column name or string.}
\item{gvar}{First treatment period, supplied as a column name or string.
Use zero for never-treated units. Must be constant within units.}
\item{cvar}{Optional one-sided formula for covariates, such as \code{~ x1 + x2}.}
\item{fevar}{Optional one-sided formula for additional fixed effects, such as
\code{~ exporter^year + importer^year}.}
\item{never}{Use never-treated controls and estimate pre-treatment effects.
Otherwise use not-yet-treated observations.}
\item{anticipation}{Positive integer giving the baseline distance in period
steps. One uses the period immediately before treatment; two moves it back one step.}
\item{hettype}{Treatment-effect specification: \code{timecohort} allows a
separate effect for each cohort-period cell; \code{time} groups effects by
period and pre/post status; \code{cohort} groups them by cohort and pre/post
status; \code{event} groups them by time relative to treatment; \code{twfe}
pools them into pre/post effects. Pre/post status is relative to the baseline
set by \code{anticipation}.}
\item{data}{A data frame, supplied as \code{data = ...}.}
\item{subset}{Logical vector aligned with \code{data}, column name or one-sided
formula selecting observations.}
\item{weights}{Positive estimation weights, treated as probability weights.
Accepts a vector aligned with \code{data}, column name or one-sided formula.}
\item{offset}{Offset in the Poisson log mean. Accepts a vector aligned with
\code{data}, column name or one-sided formula. Defaults to zero.}
\item{cluster}{Clustering formula, such as \code{~ id} or \code{~ exporter + importer}.
Defaults to the unit identifier.}
\item{...}{Additional arguments to \code{fixest::fepois()}, such as
\code{glm.tol}, \code{glm.iter} and \code{fixef.tol}. Core arguments already
set by \code{jwdid}, including \code{vcov} and \code{ssc}, must not be supplied here.}
}
\details{
Unit and time fixed effects are added automatically. Treatment persists once
it starts. Units already treated at their first observed period are removed.
If there are no never-treated units, the last cohort serves as the control and
observations at or after its treatment period are removed.

The period step is the smallest positive spacing in the usable cohort and
time support. For a cohort first treated at \code{g}, the baseline is
\code{g - anticipation * step}. With \code{anticipation > 1}, units whose first
observation is at or after \code{g - (anticipation - 1) * step} are also removed.

Observations with missing required variables and separated observations are
excluded from estimation.
The PPML path is implemented; varying-slope fixed effects are unsupported.
}
\value{
A \code{fixest} object. Use \code{summary()} for coefficient details,
\code{coef()} for coefficients, \code{vcov()} for their covariance matrix,
\code{fitted()} for fitted conditional means, and \code{estat()} for treatment
effect aggregation.
}
\seealso{\code{\link{estat}}, \code{\link[fixest]{fepois}}}
\examples{
data(gravity)
fit <- jwdid(trade, ivar = exporter^importer, tvar = year, gvar = first_year,
  fevar = ~ exporter^year + importer^year, never = TRUE, data = gravity)
att <- estat(fit, "simple")
print(att)
}
