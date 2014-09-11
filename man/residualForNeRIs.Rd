\name{residualForNeRIs}
\alias{residualForNeRIs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Residuals from prediction
%%  ~~function to do ... ~~
}
\description{
Given a model and a new data set, this function will return the residuals of the predicted values.
If the model is COX then it will return the Matingale residuals.
}
\usage{
residualForNeRIs(object, newdata,outcome,eta=0.05 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  The model used to predict the new data
%%     ~~Describe \code{object} here~~
}
  \item{newdata}{
  The new data to be predicted
%%     ~~Describe \code{newdata} here~~
}
  \item{outcome}{
 The outcome to be used to get the residuals
}
  \item{eta}{
 The for COX fitting it weights the contribution of Matingale residuals and classification residuals.
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{
Returns a vector with the differences between the prediction and the outcome
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
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Generation}
