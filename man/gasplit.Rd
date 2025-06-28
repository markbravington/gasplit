\name{gasplit}
\alias{gasplit}
\alias{gasplit2}
\alias{gasplit_missy}
\title{Flexible models for 2-category proportions from continuous data}
\description{
Suppose your data each come from one or the other of two categories (say, Eastern or Western Atlantic). You don't know for sure which one, but you do have a "measurement" for each one (e.g. some kind of log-likelihood ratio, but it could be anything) which tends to be high for one category and low for the other, with an a-priori-known distribution for each category (say, based on other samples of known category).

Now you want to estimate the proportion of category-1 samples in your data, but stratified according to covariates of each sample (say, year and location)--- and/or simply dependent on individual covariates, (say, size). Ideally, you'd like to be able to specify a flexible model formula just as you might in a GLM or GAM. Well, now you can! For example:

\preformatted{
  addmod <- gasplit( DISCRIM ~ factor(year) + factor(region) - 1,
      data= ew_data,
      d1= known_distro_of_measurement_for_category_1,
      d2= ditto_for_2
    )
}

to fit a two-way model with no interactions, to a dataset with covariates "year" and "region" based on a measurement "DISCRIM".

Having fitted your model, you can make a subsequent call to predict the proportion for (presumably) new covariate values, via the \code{predict_from_previous} argument. And, provided that your "measurement" really is a log-likelihood ratio, you can also calculate a posterior probability (of being category-1) for each case, using \code{\link{posterior}}.

If you want to use smoothers, i.e. \code{s()} terms a la \code{mgcv} (including non-sparse random effects handled via \code{s(bs="re")}), then you should call \code{gasplit2} instead. The latter requires package \pkg{RTMB}, in order to estimate smoothing parameters/random-effect variances. Results should be the same as \code{gasplit} if there are no smooth terms nor random effects; IDNK which routine is faster. NDIRC (Nor Do I Really Care.) If your covariates have missing values (MAR), you can do "full information maximum likelihood" using \code{gasplit_missy}, described next.

See \bold{Examples} of \code{\link{test_gasplit}} (qv) for practical stuff. I should update that for smoothers and missy.
\subsection{Missing covariates}{
\code{gasplit_missy} can handle Missing-At-Random (MAR) values in covariates via "full information maximum likelihood", though \emph{you} have to have already done your homework by fitting sub-models that predict the probabilities of the different possible covariate values for each case. The covariate(s) should be intrinsically discrete, or reasonably discretizable; if they are continuous, a good approximation would be use a few quantiles of their conditional distribution.

Your data will be a copy of the original dataset with extra rows added (one or more extra rows each case with missing covariates), and two extra columns, by default assumed to be called \code{MULTIMP} and \code{PROBIMP} though you can change that. Their meaning (which is very simple) is explained in the next paragraph. You might find it helpful to look at an example in \code{\link{example_prep4missy}} (qv), and the test dataset \code{\link{babydat}}. (It's not for use in actually fitting \code{gasplit_missy}--- just showing how the data prep works.)

\code{MULTIMP}--- a misnomer really, since this is actually FIML not "multiple implementation"--- should hold the "case number", so that a missing-covariate case will appear several times with the same value for some of its columns (including \code{MULTIMP} and the response variable), but different values for the missing covariate(s) in each appearance (which I might loosely call each "imputation"). The probabilities of those alternative covariate values are in \code{PROBIMP}, from that submodel(s) you've fitted (which will usually be conditional on the covariates that \emph{are} observed). Cases with no (relevant) missing covariates should have a unique \code{MULTIMP} so they only appear once, and a \code{PROBIMP} of 1. \code{MULTIMP} does not have to contain numbers; any unique case-identifier will do. Fitted values are only returned for the first occurence of each case, but take into account all "imputations" of that case. Predictions on new data don't allow for missing covariates, and the new data does not need \code{MULTIMP} or \code{PROBIMP} fields.

If there are cases where several covariates have missing values, then in principle you should produce extra rows for all combinations of possible values for that case. That could get messy for you, but fortunately it's not \emph{my} problem; \emph{you} might decide to do Multiple-Imputation instead (though it is less statistically efficient, and more fiddly).
}
\subsection{Disclaimer}{
I've tried to support most \code{mgcv} features that I more-or-less understand. For example, I \emph{think} \code{gasplit2} can cope with \code{s(...,id=...)}. However:

\itemize{
\item I'm not sure how robust my solutions are;
\item there are lots of things in \code{mgcv} that I \emph{don't} understand (and haven't used) which may not work at all.
}

An example of something that might well misfire, is trying to fix some smoothing parameters at chosen values. Let me know how y'all go...

There's deffo no support for sparsity.
}
}
\usage{
gasplit(
  formula,
  data,
  d1, d2,
  start=0.001,
  predict_from_previous= NULL,
  dbeta= FALSE
)
gasplit2(
  formula,
  data,
  d1, d2,
  predict_from_previous= NULL,
  dbeta= FALSE,
  ... # for gam()
)
gasplit_missy(
  formula,
  data,
  d1, d2,
  link_field = "MULTIMP", prob_field = "PROBIMP",
  predict_from_previous= NULL,
  dbeta= FALSE,
  ... # for gam()
)
}
\arguments{
\item{ formula}{a classic R{} model formula (eventually allowing \code{mgcv} extensions)}
\item{ data}{a dataframe...}
\item{ d1, d2}{ these are \emph{functions} giving the PDF of the response variable for category-1 and category-2 samples. You gotta pre-fit those to other data. See \bold{Examples}.}
\item{ start}{(just \code{gasplit}, not \code{gasplit2}). Either a scalar (probably close to 0) which will be used for all coefficients, or a vector with length \code{ncol(<model matrix>)}. TBH I'd leave it blank; the default will prolly work.}
\item{ link_field, prob_field}{what names the missing-data columns columns have. Probably simpler to make sure your dataset has them under the default names.}
\item{ predict_from_previous}{Optional result from a previous call to \code{gasplit}. If set, then no fitting will happen, but \code{gasplit} will make predictions at the (presumably new) values of \code{data}, based on the parameters of the previous fit.}
\item{ dbeta}{whether to calculate derivatives of each prediction WRTO the parameters (mainly for subsequent variance calculations).}
\item{ ...}{for \code{gasplit2} and \code{gasplit_missy} only, additional args to \code{mgcv::gam} (qv), eg \code{knots}.}
}
\value{
If \code{predict_from_previous} is supplied (and \code{data}, of course), all three functions return predicted proportions at each row of \code{data}. If \code{dbeta==TRUE} then the result will be not a vector, but a matrix whose first column is the point estimate, the remaining columns being its derivatives WRTO the elements of \code{beta}. You can use that to get standard errors on the predictions, something like this:
%#
\item{  t( ppn[,-1]) \%**\% V_beta \%**\% ppn[,-1]}{}
Otherwise, the functions return a list, with various elements including at least:
\item{ beta, SE_beta, V_beta}{for the raw paraameter estimates;}
\item{ ppn}{estimated proportions, casewise (or a matrix if \code{dbeta==TRUE});}
\item{ G}{setup for the underlying \code{\link{gam}}, so eg \code{<obj>$G$X} is the design matrix}
plus various convergence diagnostics. \code{gasplit2} and \code{gasplit_missy} also returns \code{$obj} which is the \code{RTMB} object (from \code{MakeADFun}) after fitting, plus the estimated "outer parameters" (smoothing parameters and/or random-effect variances) as \code{outerpars} (NB probably on log scale). The latter are treated as fixed when computing the coefficient Hessian, i.e. there's no Bayes-Empirical-Bayes correction.
\code{gasplit_missy} also returns the \code{link_field} and \code{prob_field} arguments, plus a 3-column matrix \code{$split}, with columns "PPN", <link>, and <prob> (the latter two set by \code{link_field} and \code{prob_field}). This corresponds row-wise to the input dataframe (and to \code{$G$X}), so that a single case might be split into several imputations. It is used by \code{\link{posterior}}. Note that the casewise proportions \code{$ppn} has just row per case.
}
\details{
Given parameters \code{beta} and a model-matrix \code{X} from the formula, the log-likelihood of observation \code{i} with response \code{y_i} is:

\preformatted{
  log( ppn_i*d1( y_i) + (1-ppn_i)*d2( y_i))
}

where \code{ppn_i = inv.logit( sum_j{ X_{i,j} *beta_j})}. The log-likelihood extracts the constant factor d2(y_i), to leave (in the non-missy case)

\preformatted{
  log( 1 + ppn_i*(d1_i/d2_i-1)) + const
}

In a general setting, when a covariate \code{z_i} is missing for case \code{i} but other covariates \code{x_i} may be observed, the FIML log-lik of an observation \code{y_i} is

\preformatted{
  log( sum_z{ Pr[y_i|x_i,z] * Pr[Z_i=z|x_i]})
}

where the sum ranges over possible values \code{z} that the random variable \code{Z_i} might take, and the second prob is your "submodel". The point here is that the weighted sum is (and should be) applied to the \emph{unlogged} probabilities, \emph{before} the whole lot gets logged finally. That's why the code needs to be adapted.

For our particular application, this becomes:

\preformatted{
  log( sum_z{ (ppn_i(z)*d1( y_i) + (1-ppn_i(z))*d2( y_i)) * prob(z_i=z)}
  = log( 1 + (d1_i/d2_i-1) * sum_z{ ppn_i(z) * prob( z_i=z)}
}

with some obvious & harmless abuse of notation.

Smoothers and random-effects add "penalties" (ie priors) on the coefficients in the usual \code{mgcv} way, and \code{RTMB} computes the automated Laplace approximation required for consistent estimation of smoopars/RE-variances.
}
\examples{
?test_gasplit
?example_prep4gasplit
}
\seealso{test_gasplit, example_prep4gasplit
}
\keyword{misc}
