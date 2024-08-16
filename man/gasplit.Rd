\name{gasplit}
\alias{gasplit}
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

See EXAMPLES of \code{\link{test_gasplit}} (qv) for practical stuff.
\subsection{To do}{
Move this to RTMB, and add random effects, so that you could write eg

\preformatted{
  fancymod <- gasplit( LGLR ~ year + region + s( year, region, bs='re'), ...)
}

using \code{mgcv} syntax to put a random-effect on each interaction term. It's easy enough in principle, but everything takes time...
}
}
\usage{
gasplit(formula, data, d1, d2, start=0.001,
    predict_from_previous=NULL)
}
\arguments{
\item{ formula}{a classic R{} model formula (eventually allowing \code{mgcv} extensions)}
\item{ data}{a dataframe...}
\item{ d1, d2}{ these are \emph{functions} giving the PDF of the response variable for category-1 and category-2 samples. You gotta pre-fit those to other data. See EXAMPLES.}
\item{ start}{Either a scalar (probably close to 0) which will be used for all coefficients, or a vector with length \code{ncol(<model matrix>)}}
\item{ predict_from_previous}{Optional result from a previous call to \code{gasplit}. If set, then no fitting will happen, but \code{gasplit} will make predictions at the (presumably new) values of \code{data}, based on the parameters of the previous fit.}
}
\value{A list with elements \code{par} (raw parameter estimates from \code{optim}) and \code{ppn}, which is a vector of estimated proportion-1 for each \emph{observation}. If the "strata" are discrete, there will be lots of identical \code{ppn} values, but not necessarily the case.
}
\details{
Given parameters \code{beta} and a model-matrix \code{X} from the formula, the log-likelihood of observation \code{i} with response \code{y_i} is:

\preformatted{
  ppn_i*d1( y_i) + (1-ppn_i)*d2( y_i)
}

where \code{ppn_i = inv.logit( sum_j{ X_{i,j} *beta_j})}.
}
\keyword{misc}
