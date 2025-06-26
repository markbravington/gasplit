\name{posterior}
\alias{posterior}
\title{Posterior probability}
\description{
After fitting a \code{\link{gasplit}} model to estimate the \emph{prior} probability of whateverness by covariates, you can also calculate a per-observation \emph{posterior} probability.

The posteriors are also estimates, since they depend on the estimated model coefficients \code{beta}. You can also get their derivatives WRTO the coefficients (useful for downstream variance propagation).

Missing data can be handled a la \code{\link{gasplit_missy}} (qv) by setting the \code{link_field} and \code{prob_field} arguments. This means the outputs are aggregated to case level with weighted sums over the imputations (which I think is the correct approach). It is most likely to be useful for extracting posteriors from an existing model that was fitted by \code{\link{gasplit_missy}} (when it will be automatic), rather than for brand-new predictions. But, whatever floats your modeling boat.
}
\usage{
#posterior(
#  object,
#  newdata = NULL,
#  dbeta = FALSE.
#  link_field= object$link_field,
#  prob_field= object$prob_field
#)
}
\arguments{
\item{ object}{from previous \code{gasplit} etc call}
\item{ newdata}{default is to use the original data, but it can accept a dataframe with the necessary covariates. If \code{newdata} includes the columns named in \code{link_field} and \code{prob_field}, then results will be aggregated casewise as per \code{gasplit_missy}.}
\item{ dbeta}{if TRUE, return derivs of the posterior probs WRTO the params}
\item{ link_field, prob_field}{optional for "multiple imputation" with missing data; see \code{gasplit_missy}. Unlikely with de novo \code{newdata} but you never know.}
}
\value{If \code{dbeta==FALSE}, a vector of estimated posterior probs, one per case. If not, then a matrix where the first column is the posteriors, and the remaining columns are the derivs WRTO each component of \code{beta}.
}
\seealso{\code{\link{gasplit}}, \code{\link{gasplit_missy}}
}
\examples{
sim2 <- make_fake_ppns(
    meanE= +2, # not-brilliant separation
    mean_nsamp=200, # reasonable samp sizes
    interpow=0) # no interactoins
hist( sim2$LGLR, nc=40)
test <- test_gasplit( sim2)
with( test, plot( tru_ppnE, fit_ppnE))
abline( 0,1)
posty <- posterior( test$gg1, dbeta=TRUE)
plot( sim2$LGLR, logit( posty[,1]), pch='.')
abline( 0,1,col='red')
# ... Points with high LGLR tend to pushed upwards in the posterior
# ... which makes sense, coz if you've got high LGLR, your neighbours prolly will too
# ... so your prior would be pushed up
}
\keyword{misc}
