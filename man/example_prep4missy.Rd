\name{example_prep4missy}
\alias{example_prep4missy}
\alias{babydat}
\title{Schematic data-prep for 'gasplit_missy'}
\description{
\code{example_prep4missy} is a "template" showing how to prep a dataset with (one) sometimes-missing covariate (ie field) into a form suitable for \code{\link{gasplit_missy}} (qv). You could adapt \code{example_prep4missy} to your own situation.

There is an example dataset \code{babydat} just to demonstrate the prep (you can't use it to fit a model, coz the \code{d1/d2} distributions are missing). The covariate in question is \code{Sex}, taking values \code{"F"}, \code{"M"}, or \code{NA}. Someone has already fitted a model which predicts the probability that a case will be Female (\code{"F"}) based on its other covariates (but \bold{not} on whatever the ultimate response variable will be!); that prediction has been stored in the field \code{PrFem_yl}. Cases where \code{Sex} is known are left untouched; cases where \code{Sex} is \code{NA} are duplicated (since there are only two possible values for this covariate), with one copy getting \code{Sex="F"}, and the other \code{Sex="M"}; the field \code{PROBIMP} is set accordingly based on \code{PrFem_yl}. The \code{MULTIMP} field is used to link together the multiple versions of each case.

Remember that \emph{you} have to develop whatever submodel for probabilities of the alternative values that \emph{your} missing covariates might take! If you've got missing values on continuous covariate, rather than a discrete one like Sex, my recommendation is to add rows based on a few quantiles of the covariate (which will potentially be at different values for each case), and give them all equal weight (ie in the analog of \code{PrFem_yl}). I'd start by trying 5 quantiles (0.1; 0.3; ...; 0.9), but you could increase & see if it matters.
}
\usage{example_prep4missy( data)
}
\arguments{\item{ data}{data.frame as per description.}
}
\seealso{\code{\link{gasplit_missy}}
}
\examples{
babydat      # with some missing values for 'Sex'
example_prep4missy( babydat)
}
\keyword{misc}
