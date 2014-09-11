\name{improvedResiduals}
\alias{improvedResiduals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates the significance of the reduction of predicted residuals
%%  ~~function to do ... ~~
}
\description{
Given a set of two residuals (new vs. old), this function will test the hypothesis that the new residuals are better than the old residuals.
Four p-values are provided: Paired t-test, Wilcoxon signed rank test, The binomial sing-test and the F test
The proportion of subjects that improved their residuals, the proportion that worsen their residuals and the net residual improvement will be returned
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
improvedResiduals(oldResiduals,newResiduals,
testType=c("Binomial","Wilcox","tStudent","Ftest"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oldResiduals}{
  the residuals of the original model
}
  \item{newResiduals}{
  the residuals of the new model
}
  \item{testType}{
  The test to be used to provide the final p.value.
}

}
\details{
This function can be used to test the hypothesis that the new residuals are "better" than the old residuals.
Four type of tests are provided.
The paired t-test compares the absolute value of the residuals. 
The paired wilcox test compares the absolute value of residuals
The sign test test if the number of subjects that improve residuals is greater than the number of subjects that worsen residuals
The F test is the standard test for testing that the residual variance is better in the new residuals.
The proportion of subjects that improved residuals,and the ones that worst residuals is returned.
The net proportion of improvement in residuals is returned.
%%  ~~ If necessary, more details than the description above ~~
}
\value{

%%  If it is a LIST, use
  \item{p1}{Proportion of subjects that improved their residuals}
  \item{p2}{Proportion of subjects that worsen their residuals}
  \item{NeRI}{The net improvement (p1-p2)}
  \item{binom.pValue}{The p value associated with a significant improvement in residuals}
  \item{wilcox.pValue}{The single sided p value of the wilcox test comparing the absolute value of the new vs. old residuals}
  \item{t.test.pValue}{The single sided p value of the t-test comparing the absolute value of the new vs. old residuals}
  \item{F.test.pValue}{The single sided p value of the F-test comparing the residual variances}
  \item{p.Value}{The one tail p value of the specified test}
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
%%
%%
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Generation}
