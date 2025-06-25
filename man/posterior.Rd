\name{posterior}
\alias{posterior}
\title{Posterior probability}
\description{
After fitting a \code{\link{gasplit}} model to estimate the \emph{prior} probability of whateverness by covariates, you can calculate a per-observation \emph{posterior} probability.

The posteriors are also estimates, since they depend on the estimated model coefficients \code{beta}. You can also get their derivatives WRTO the coefficients (useful for downstream variance propagation).

There's no attempt to predict "properly" for multiple imputations. You can DIY by predicting for each case, and manually weighted-summing the predictions; the mean will be OK. However, uncertainty won't be fully propagated. I reckon that's a NFN...
\subsection{Requirements}{\code{posterior} \bold{assumes} that your discrimination statistic is actually a log-likelihood ratio--- not just some arbitrary stat that you made up, not a LGLR multiplied by 3, etc... and \bold{WATCH OUT} that the sign is correct! (ie that it's lglk(catgeory-1 / cat-2)) not the other way round... be careful!!! (Check that there's a \bold{positive} correlation between prior and posterior, and between lglk and posterior...)
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
