\name{summaryReport}
\alias{summaryReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reports the Univariate analysis, the cross-validation analysis and the correlation analysis .
%%  ~~function to do ... ~~
}
\description{
The function takes the variables of the cross-validation analysis and extracts its univariate analysis and its correlation analysis.
Then it prints the cross-validation, the univariate analysis and the correlated variables.
At the output returns a list of each one.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

summaryReport(univariateObject,
summaryBootstrap,
listOfCorrelatedVariables=NULL,
digits=2
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{univariateObject}{
  A data frame that contains the results of the univariate analysis

}
  \item{summaryBootstrap}{
  An summary.bootstrappedIDI object. It will extract the performance table and the coefficient table from the summary
}
  \item{listOfCorrelatedVariables}{
  A data frame that contains the list of significantly correlated variables for each variable
}
  \item{digits}{
  The number of digits to be used in the print function
}

}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{

%%  If it is a LIST, use
  \item{performance.table}{A matrix data frame that list the AUC, accuracy, sensitivity, specificity of the bootstrapped model }
  \item{coefStats}{A matrix that list all the model features along with its univariate statistics and bootstrapped coefficients}
  \item{cor.varibles}{List the features that are correlated to the model variables}
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
\keyword{Model_Inspection}
