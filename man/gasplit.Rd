\name{gasplit}
\alias{gasplit}
\alias{gasplit2}
\alias{gasplit_multimp}
\title{Flexible models for 2-category proportions from continuous data}
\description{
Suppose your data each come from one or the other of two categories (say, Eastern or Western Atlantic). You don't know for sure which one, but you do have a "measurement" for each one (say, a log-likelihood ratio) which tends to be high for one category and low for the other, with an a-priori-known distribution for each category (say, based on other samples of known category).

Now you want to estimate the proportion of category-1 samples in your data, but stratified according to covariates of each sample (say, year and location)--- and/or simply dependent on individual covariates, (say, size). Ideally, you'd like to be able to specify a flexible model formula just as you might in a GLM or GAM. Well, now you can! For example:

\preformatted{
  addmod <- gasplit( LGLR ~ factor(year) + factor(region) - 1,
      data= ew_data,
      d1= known_distro_of_measurement_for_category_1,
      d2= ditto_for_2
    )
}

to fit a two-way model with no interactions, to a dataset with covariates "year" and "region" based on a measurement "LGLR".

Having fitted your model, you can make a subsequent call to predict the proportion for (presumably) new covariate values, via the \code{predict_from_previous} argument..

If you want to use smoothers, i.e. \code{s()} terms a la \code{mgcv} (including non-sparse random effects handled via \code{s(bs="re")}), then you should call \code{gasplit2} instead. The latter requires package \pkg{RTMB}, in order to estimate smoothing parameters/random-effect variances. Results should be the same as \code{gasplit} if there are no smooth terms nor random effects; IDNK which routine is faster. NDIRC (Nor Do I Really Care.) If your covariates have missing values and you want multiple imputation, see next subsection for \code{gasplit_multimp}.

See \bold{Examples} of \code{\link{test_gasplit}} (qv) for practical stuff. I should update that for smoothers and multimp.
\subsection{Multiple imputation}{\code{gasplit_multimp} can handle missing values in discrete covariates, though \emph{you} have to have already done your homework by fitting sub-models that predict the probabilities of the different possible covariate values for each case. Your data will need two extra columns, by default assumed to be called \code{MULTIMP} and \code{PROBIMP} though you can change that. \code{MULTIMP} holds the "case number", so that a multiply-imputed case will appear several times with the same value for some of its columns (including \code{MULTIMP} and the response variable), but different values for the missing covariate(s) in each imputation. The probabilities of those alternative imputed values are in \code{PROBIMP}, from that submodel(s) you've fitted. Cases with no (relevant) missing covariates should have a unique \code{MULTIMP} so they only appear once, and a \code{PROBIMP} of 1. \code{MULTIMP} does not have to contain numbers; any unique case-identifier will do. Fitted values are only returned for the first occurence of each case, but are based on all imputations of that case. Predictions on new data don't allow for missing covariates, and the new data does not need \code{MULTIMP} or \code{PROBIMP} fields.
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
gasplit(formula, data, d1, d2, start=0.001,
    predict_from_previous=NULL)
gasplit2(
  formula,
  data,
  d1, d2,
  predict_from_previous= NULL,
  ... # for gam()
)
gasplit_multimp(
  formula,
  data,
  d1, d2,
  link_field = "MULTIMP", prob_field = "PROBIMP",
  predict_from_previous= NULL,
  ... # for gam()
)
}
\arguments{
\item{ formula}{a classic R{} model formula (eventually allowing \code{mgcv} extensions)}
\item{ data}{a dataframe...}
\item{ d1, d2}{ these are \emph{functions} giving the PDF of the response variable for category-1 and category-2 samples. You gotta pre-fit those to other data. See \bold{Examples}.}
\item{ start}{(just \code{gasplit}, not \code{gasplit2}). Either a scalar (probably close to 0) which will be used for all coefficients, or a vector with length \code{ncol(<model matrix>)}. TBH I'd leave it blank; the default will prolly work.}
\item{ link_field, prob_field}{what names the multiple-imputation columns have.}
\item{ predict_from_previous}{Optional result from a previous call to \code{gasplit}. If set, then no fitting will happen, but \code{gasplit} will make predictions at the (presumably new) values of \code{data}, based on the parameters of the previous fit.}
\item{ ...}{for \code{gasplit2} and \code{gasplit_multimp} only, additional args to \code{mgcv::gam} (qv), eg \code{knots}.}
}
\value{
If \code{predict_from_previous} is supplied (and \code{data}, of course), all three functions return predicted proportions at each row of \code{data} (without SEs--- I prolly orta add that).
Otherwise, the functions return a list, with various elements:
\item{ beta, SE_beta, V_beta}{for the raw paraameter estimates;}
\item{ ppn}{estimated proportions, casewise;}
\item{ G}{setup for the underlying \code{\link{gam}}, so eg \code{G$X} is the design matrix}
plus various convergence diagnostics. \code{gasplit2} and \code{gasplit_multimp} also returns \code{obj} which is the \code{RTMB} object (from \code{MakeADFun}) after fitting, plus the estimated "outer parameters" (smoothing parameters and/or random-effect variances) as \code{outerpars} (NB probably on log scale). The latter are treated as fixed when computing the coefficient Hessian, i.e. there's no Bayes-Empirical-Bayes correction.
}
\details{
Given parameters \code{beta} and a model-matrix \code{X} from the formula, the log-likelihood of observation \code{i} with response \code{y_i} is:

\preformatted{
  log( ppn_i*d1( y_i) + (1-ppn_i)*d2( y_i))
}

where \code{ppn_i = inv.logit( sum_j{ X_{i,j} *beta_j})}. The log-likelihood extracts the constant factor d2(y_i), to leave (in the non-multimp case)

\preformatted{
  log( 1 + ppn_i*(d1_i/d2_i-1)) + const
}

For multiple imputation of one or more covariates \code{z}, this becomes

\preformatted{
  log( sum_z{ (ppn_i(z)*d1( y_i) + (1-ppn_i(z))*d2( y_i)) * prob(z_i=z)}
  = log( 1 + (d1_i/d2_i-1) * sum_z{ ppn_i(z) * prob( z_i=z)}
}

with some obvious & harmless abuse of notation.

Smoothers and random-effects add "penalties" (ie priors) on the coefficients in the usual \code{mgcv} way, and \code{RTMB} computes the automated Laplace approximation required for consistent estimation of smoopars/RE-variances.
}
\examples{
?test_gasplit
}
\seealso{test.gasplit
}
\keyword{misc}
