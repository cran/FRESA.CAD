\name{summary.bootstrapValidation}
\alias{summary}
\alias{summary.bootstrapValidation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a report of the result of the bootstrap-based validation: bootstrapValidation
%%  ~~function to do ... ~~
}
\description{
The summary function prints the results of the bootstrap-based-validation of binary classification models
The first table reports the accuracy, sensitivity, specificity and AUC for the training set with is confidence intervals.
The test accuracy, sensitivity, specificity and AUC are also reported.
The second table reports the model coefficients and its corresponding IDI, NRI performance.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
\method{summary}{bootstrapValidation}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A bootstrapValidation object
}
  \item{...}{ parameters to print
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{

%%  If it is a LIST, use
  \item{performance}{A data matrix with the results of the bootstrapping procedure}
  \item{summary}{A summary of the model as reported by the native fitting function}
  \item{coef}{A data matrix with the coefficients, their 95 bootstrapped CI, 
  the IDI, the 95 CI of IDI, the NRI, the 95 CI of NRI}
  \item{performance.table}{A data matrix with the tabulated results of the blind/test accuracy, sensitivity, specificities and ROC AUC}  
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
\keyword{Model_Inspection}
