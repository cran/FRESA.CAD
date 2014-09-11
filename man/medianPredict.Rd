\name{medianPredict}
\alias{medianPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The median prediction from a list of models
%%  ~~function to do ... ~~
}
\description{
Given a list of models formulas, it will train them and it will return the median prediction on a new data set
It also provides the KNN prediction using the features listed in the models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
medianPredict(formulaList,traindata,newdata=NULL, 
predictType=c("prob","linear"),
modelType = c("LOGIT", "LM","COX"),Outcome="CLASS",k=0,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formulaList}{
   The list of formulas to be fitted and predicted
}
  \item{traindata}{
  The data to train the models
}
  \item{newdata}{
  The data set to be predicted
}
  \item{predictType}{
   the type of prediction: probability or linear
}
  \item{modelType}{
   the type of fitting: LOGIT, LM or COX
}
  \item{...}{
   additional parameters for fitting a default glm model
}
  \item{Outcome}{
   the outcome for KNN
}
  \item{k}{
   the number of neighbours for KNN. If KNN is less than 0 it will not perform the KNN
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{

%%  If it is a LIST, use
  \item{modelPredict}{The median of the predicted outcome}
  \item{KNN.Predict}{The median of the KNN predicted outcome}
  \item{predictions}{A data matrix with the prediction of each model}
  \item{KNNpredictions}{A data matrix with the prediction of each KNN}
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
\keyword{Model_Inspection}
