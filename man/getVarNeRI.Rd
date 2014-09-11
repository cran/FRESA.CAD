\name{getVarNeRI}
\alias{getVarNeRI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyses each model term by comparing the residuals of the full vs the reduced model.
%%  ~~function to do ... ~~
}

\description{
Provides an analysis of the effect of each model term by comparing the residuals of the full model vs the reduced model.
The model is fitted using the train data set, but analysis of residual improvement is done on the training and test sets.
Full and reduced-model residuals are compared by paired t-test, paired Wilcox test, binomial sign test and the F test on residual variance.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
getVarNeRI(object, dataframe, Outcome = "Class", 
type = c("LM","LOGIT", "COX"),testdata=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  The model to be tested
%%     ~~Describe \code{object} here~~
}
  \item{dataframe}{
  The train data frame to be used to fit the model

}
  \item{Outcome}{
  The name of the outcome column to be predicted
%%     ~~Describe \code{Outcome} here~~
}
  \item{type}{
  The type of fit : COX, LOGIT or LM

}
  \item{testdata}{
  An additional test data frame for reporting test predictions

}

}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{

%%  If it is a LIST, use
  \item{tP.value}{The p value of the t-test}
  \item{BinP.value}{The p value of the binomial test}
  \item{WilcoxP.value}{The p value of the Wilcoxon test}
  \item{FP.value}{The p value of the F test}
  \item{NeRIs}{The proportion of subjects that improved minus the subjects that worsen residuals}
  \item{testData.tP.value}{The p value of the t-test on the test data}
  \item{testData.BinP.value}{The p value of the binomial test on the test data}
  \item{testData.WilcoxP.value}{The p value of the Wilcoxon test on the test data }
  \item{testData.FP.value}{The p value of the F test on the test data }
  \item{testData.NeRIs}{The proportion of subjects that improved minus the subjects that worsen residuals on the test data}
%% ...
}

%%\references{
%% ~put references to the literature/web site here ~
%%}

\author{
Jose G. Tamez-Pena

}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

%%\examples{
%%
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Inspection}
