\name{predictForFresa}
\alias{predictForFresa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear or Probabilistic prediction from models
%%  ~~function to do ... ~~
}
\description{
This predict function returns the observed predicted values of the candidates models.
The models generate the linear prediction and the probabilistic values are generated using the logistic transform
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predictForFresa(object, newdata, type = c("prob", "linear"))
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
  \item{type}{
  The type of prediction: Linear or prob.
  The prob type computes the logisitc transformation of the linear predictor

}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{
Returns a vector with the predicted values
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
