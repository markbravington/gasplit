\name{gasplit}
\alias{gasplit}
\alias{gasplit2}
\title{Flexible models for 2-category proportions from continuous data}
\description{
Suppose your data each come from one or the other of two categories (say, Eastern or Western Atlantic). You don't know for sure which one, but you do have a "measurement" for each one (say, a log-likelihood ratio) which tends to be high for one category and low for the other, with an a-priori-known distribution for each category (say, based on other samples of known category).

Now you want to estimate the proportion of category-1 data in your data, but \emph{stratified} according to covariates of the data (say, year and location)--- and/or simply dependent on individual covariates, (say, size). Ideally, you'd like to be able to specify a flexible model formula just as you might in a GLM or GAM. Well, now you can! For example:

\preformatted{
  addmod <- gasplit( LGLR ~ year + region - 1, data= ew_data,
      d1= known_distro_of_measurement_for_category_1,
      d2= ditto_for_2
    )
}

to fit a two-way model with no interactions, to a dataset with covariates "year" and "region" based on a measurement "LGLR".

Having fitted your model, you can make a subsequent call to predict the proportion for (presumably) new covariate values, via the \code{predict_from_previous} argument..

If you want to use smoothers, i.e. \code{s()} terms a la \code{mgcv}, then you should call \code{gasplit2} instead. The latter requires \code{RTMB}. Results should be the same as \code{gasplit} if there are no smooth terms; IDNK which one is faster. NDIRC (Nor Do I Really Care.)

See \bold{Examples} of \code{\link{test_gasplit}} (qv) for practical stuff.
}
\usage{
gasplit(formula, data, d1, d2, start=0.001,
    predict_from_previous=NULL)
gasplit2(
  formula,
  data,
  d1, d2,
  predict_from_previous= NULL,
  ... # for gam()
)
}
\arguments{
\item{ formula}{a classic R{} model formula (eventually allowing \code{mgcv} extensions)}
\item{ data}{a dataframe...}
\item{ d1, d2}{ these are \emph{functions} giving the PDF of the response variable for category-1 and category-2 samples. You gotta pre-fit those to other data. See \bold{Examples}.}
\item{ start}{(just \code{gasplit}, not \code{gasplit2}). Either a scalar (probably close to 0) which will be used for all coefficients, or a vector with length \code{ncol(<model matrix>)}. TBH I'd leave it blank; the default will prolly work.}
\item{ predict_from_previous}{Optional result from a previous call to \code{gasplit}. If set, then no fitting will happen, but \code{gasplit} will make predictions at the (presumably new) values of \code{data}, based on the parameters of the previous fit.}
\item{ ...}{for \code{gasplit2} only, additional args to \code{mgcv::gam} (qv), eg \code{knots}.}
}
\value{
If \code{predict_from_previous} is supplied (and \code{data}, of course), \code{gasplit} and \code{gasplit2} both return predicted proportions at each row of \code{data}.
Otherwise, \code{gasplit} returns a list with elements \code{par} (raw parameter estimates from \code{optim}) and \code{ppn}, which is a vector of estimated proportion-1 for each \emph{observation}. If the "strata" are discrete, there will be lots of identical \code{ppn} values, but not necessarily the case.
\code{gasplit2} returns the \code{RTMB} object (from \code{MakeADFun}) after fitting. A few extra components are added: \code{G} (from gam setup), \code{beta} (underlying params; could perhaps be used for variance computations via delta-method and hessian and blah-blah) and \code{ppn} (fitted values), plus \code{convergence}, \code{message}, and \code{evaluations} from \code{nlminb}. Calling \code{<obj>$report()} on the output will give a list with \code{beta} and \code{ppn}. If \code{gasplit2} is called without any \code{s()} terms, then the internal details of the \code{RTMB} object will be a bit different, but the extra components will be the same.
}
\details{
Given parameters \code{beta} and a model-matrix \code{X} from the formula, the log-likelihood of observation \code{i} with response \code{y_i} is:

\preformatted{
  ppn_i*d1( y_i) + (1-ppn_i)*d2( y_i)
}

where \code{ppn_i = inv.logit( sum_j{ X_{i,j} *beta_j})}.
}
\examples{
?test_gasplit
}
\seealso{test.gasplit
}
\keyword{misc}
