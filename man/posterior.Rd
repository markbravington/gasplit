\name{posterior}
\alias{posterior}
\title{Posterior probability}
\description{
After fitting a \code{\link{gasplit}} (or \code{\link{gasplit2}} or \code{\link{gasplit_missy}}) model to estimate the \emph{prior} probability of whateverness by covariates, you can calculate a per-observation \emph{posterior} probability.

The posteriors are also estimates, since they depend on the estimated model coefficients \code{beta}. You can also get their derivatives WRTO the coefficients (useful for downstream variance propagation).

If the original model used multiple imputation (\code{\link{gasplit_missy}}), then the posterior for a case is the weighted sum of posteriors for each of its imputations, and likewise for the derivs. The result is inserted only for the first imputation; ; the rest are all set to NA.
\subsection{Requirement}{\code{posterior} \bold{assumes} that your discrimination statistic is actually a log-likelihood ratio--- not just some arbitrary stat that you made up, not a LGLR multiplied by 3, etc...
}
}
\usage{posterior(object, newdata = NULL, dbeta = FALSE)
}
\arguments{
\item{ object}{from previous \code{gasplit} etc call}
\item{ newdata}{default is to use the original data, but it can al}
\item{ dbeta}{if TRUE, return derivs of the posterior probs WRTO the params}
}
\value{If \code{dbeta==FALSE}, a vector of estimated posterior probs, one per case. If not, then a matrix where the first column is the posteriors, and the remaining columns are the derivs WRTO each component of \code{beta}.
}
\seealso{\code{\link{gasplit}}, \code{\link{gasplit_missy}}
}
\examples{
# ?test_gasplit
sim2 <- make_fake_ppns(
    meanE= -2, # not-brilliant separation
    mean_nsamp=200, # reasonable samp sizes
    interpow=0) # no interactoins
hist( sim2$LGLR, nc=40)
test <- test_gasplit( sim2)
with( test, plot( tru_ppnE, fit_ppnE))
abline( 0,1)
posty <- posterior( test$gg1, dbeta=TRUE)
plot( inv.logit( sim2$LGLR), posty[,1])
}
\keyword{misc}
