\name{modelFitting}
\alias{modelFitting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits a model to the data
%%  ~~function to do ... ~~
}
\description{
A FRESA.CAD interface to fit a model to the dataframe. Fitting: Linear, Logistic or COX
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
modelFitting(model,dataframe,type = c("LOGIT", "LM","COX"),...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
   A model formula to be fitted
}
  \item{dataframe}{
  a data frame with the feature data
}
  \item{type}{
   the type of fitting: LOGIT, LM or COX
}
  \item{...}{
   additional parameters for fitting a default glm model
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

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
\value{
Returns the fitted model to the prescribed type
}

%%\examples{
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Generation}
