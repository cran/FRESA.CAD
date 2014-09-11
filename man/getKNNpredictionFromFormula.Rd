\name{getKNNpredictionFromFormula}
\alias{getKNNpredictionFromFormula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get the KNN predicted classification using the terms of a given formula
%%  ~~function to do ... ~~
}
\description{
Given a formula and a train set. This function will return the classification of the samples of the test set using a KNN algorithm with euclidean distance.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getKNNpredictionFromFormula(modelFormula,trainData,testData,outcome="CLASS",k=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelFormula}{
   The formula from witch to extract the variables
}
  \item{trainData}{
  The data to train the models
}
  \item{testData}{
  The data set to be predicted
}
  \item{outcome}{
   what to be predicted
}
  \item{k}{
   the number of neighbours
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{

%%  If it is a LIST, use
  \item{prediction}{The results of the prediction}
  \item{prob}{The proportion of K neighbours that the class}
  \item{binProb}{The proportion of K neighbours that the class = 1}
  \item{featureList}{The features used by the KNN procedure}
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

\seealso{
See Also \code{knn} from the "class" package
}

%%\examples{
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Inspection}
