\name{bootstrapValidationNeRI}
\alias{bootstrapValidationNeRI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap validation of regression models.
%%  ~~function to do ... ~~
}
\description{
bootstrapValidationNeRI bootstrap the model N times to estimate for each variable the empirical bootstrapped distribution of coefficients.
At each bootstrap the non-observed data is predicted by the trained model and statistic of test-prediction are stored and reported. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapValidationNeRI(fraction = 1, loops = 200, model.formula, 
Outcome, dataframe, type = c("LM", "LOGIT", "COX"),plots=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fraction}{
  The fraction of data to be sampled with replacement (bootstrap sampling)
}
  \item{loops}{
  The number of bootstrapped steps
}
  \item{model.formula}{
  The formula to be bootstrapped and validated. It has to be either a COX object or a LM model.
%%     ~~Describe \code{model.formula} here~~
}
  \item{Outcome}{
  The outcome variable to which the accuracy, sensitivity and specificity will be computed.
}
  \item{dataframe}{
  The training data frame 
}
  \item{type}{
  The type of fit: LM", "LOGIT", "COX"
}
  \item{plots}{
  Set to FALSE to avoid showing the density plots 
}
}
\details{
The bootstrap validation will estimate the CI of the model coefficients and the NeRI.
It will compute the test RMSE and the train RMSE as well as the distribution of the NeRI pvalues
%%  ~~ If necessary, more details than the description above ~~
}
\value{

%%  If it is a LIST, use
  \item{outcome}{A vector of the outcome}
  \item{coefficients}{The matrix of the sampled coefficients}
  \item{NeRis}{The matrix of each coefficient NeRI}
  \item{tStudent.pvalues}{the t-student p value of the NeRI}
  \item{wilcox.pvalues}{the wilcox p value of the NeRI}
  \item{bin.pvlaues}{the binomial p value of the NeRI}
  \item{F.pvlaues}{the F p value of the NeRI}
  \item{test.tStudent.pvalues}{the t-student p value of the test NeRI}
  \item{test.wilcox.pvalues}{the wilcox p value of the test NeRI}
  \item{test.bin.pvlaues}{the binomial p value of the test NeRI}
  \item{test.F.pvlaues}{the F p value of the test NeRI}
  \item{testPrediction}{the predicted observations on the test set}
  \item{testOutcome}{the true outcome of the test set}
  \item{testResiduals}{The test residuals}
  \item{trainPrediction}{The training predictions}
  \item{trainOutcome}{The training outcomes}
  \item{trainResiduals}{The training residuals}
  \item{testRMSE}{The global test RMSE}
  \item{trainRMSE}{The global train RMSE}
  \item{trainSampleRSME}{The bootstrap train sample RMSE}
  \item{testSampledRSME}{The bootstrap test sample RMSE}
%% ...
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Jose G. Tamez-Pena

}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Diagnosis}
