\name{make_fake_ppns}
\alias{make_fake_ppns}
\alias{test_gasplit}
\title{Simulate dataset for 'gasplit'}
\description{
\code{make_fake_ppns} simulates a dataset for fitting with \code{\link{gasplit}}, using two factors \code{Y} and \code{Z} and giving a response variable \code{LGLR} (so the strata are all combinations of \code{Y} and \code{Z}). You can control the number of factor-levels, distributions of \code{LGLR} for category-1 and category-2 observations (and thus the "power" of \code{LGLR} for assignment), the number of samples per stratum, the average strength of each Y-effect and Z-effect, the average strength of Y*Z interactions, etc. It requires the \pkg{offarray} package (purely for my programming convenience; it's not apparent in the output).

\code{test_gasplit} fits a \code{\link{gasplit}} model (with no interactions) to the output of \code{make_fake_ppns}, running the latter itself if required. The output is the (true and) fitted proportions per strata, which can be compared: see \bold{Examples}. Note that it's a bit harsh to simulate with \code{interpow>0} and then test with \code{\link{gasplit}}, coz the latter assumes there's no interactions. (Of course, you can fit a model with interactions if you want.)
}
\usage{
make_fake_ppns(nyears = 4, nzones = 3, interpow = 0.1,
    df_t = 5, mean_nsamp = 100,  meanE = 1, prange = 1, seed = 2)
test_gasplit( sim=NULL, ...)
}
\arguments{
\item{ nyears, nzones}{number of factor levels}
\item{ interpow}{average strength of an interaction relative to a main effect}
\item{ mean_nsamp}{actually there are \emph{exactly} this many total samples from each stratum}
\item{ df_t, meanE}{\code{df_t} degrees-of-freedom for the t-distributions of response by category, which will have variance 1 and means \code{+meanE} and \code{-meanE}. A \code{meanE} of 1 corresponds to pretty poor separation; 3 is pretty good.}
\item{ prange}{how strong the effects should be. If they are too strong, then}
\item{ seed}{random number seed. The system seed is restored on completion.}
\item{ sim, ...}{(test_gasplit) \code{sim} should come from a previous call to \code{make_fake_ppns}. If not supplied, then \code{make_fake_ppns} will be run first, using the \code{...} arguments.}
}
\value{\code{make_fake_ppns} returns a dataframe that has all the stuff you need, plus an attribute \code{truth} which is a list containing the real dope, as stratum-specific \code{offarray}s. Note that \code{truth$samp_ppnE} is the actual ppn of category-1 ("E") samples in each stratum, which will differ from the nominal value given by the factor strengths and interactions because of binomial sampling variability.
}
\seealso{\code{\link{gasplit}}
}
\examples{
if( require( offarray))\{
sim2 <- make_fake_ppns(
    meanE= -2, # not-brilliant separation
    mean_nsamp=200, # reasonable samp sizes
    interpow=0) # no interactoins
hist( sim2$LGLR, nc=40) # to see how good the separation is
test2 <- test_gasplit( sim2)
with( test2, plot( tru_ppnE, fit_ppnE))
abline( 0,1)
# Make up some new data...
splodge <- with( test2, simpure[ 1:7,])
# ... and predict. The result of the original call to gasplit() is
# stored in test2$gg1
splodge_pred <- gasplit( data=splodge, predict_from_previous=test2$gg1)
# ... for normal non-test situations, do something more like this:
if( FALSE)\{
  pfit <- gasplit( myform, mydata)
  ppred <- gasplit( data=newmydata, predict_from_previous=pfit)
\}
\}
}
\keyword{misc}
