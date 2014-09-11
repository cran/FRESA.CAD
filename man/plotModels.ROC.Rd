\name{plotModels.ROC}
\alias{plotModels.ROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the test ROC curves of each one of the Cross validated models
%%  ~~function to do ... ~~
}
\description{
Plot each one of the models found in the cross validation process.
It will aggregate the models into a single prediction performance, plotting the resulting ROC curve (model coherence ROC)
Finally it will plot the mean sensitivity for the given set of specificities.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotModels.ROC(modelPredictions, number.of.models=0,
specificities=c(0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelPredictions}{
   A data frame created by the model cross validation process
}
  \item{number.of.models}{
  the maximum number of models to plot
}
  \item{specificities}{
   list of specificities for ROC sensitivities 
}
  \item{...}{
   additional parameters for the ROC function
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

%%  If it is a LIST, use
  \item{ROC.AUCs}{A list with the individual AUC of the ROC}
  \item{mean.sensitivities}{A vector of the mean sensitivity given the specificities}
  \item{model.sensitivities}{A matrix where each row is the sensitivities for each model}
  \item{specificities}{The specificities used to calculate the sensitivities}
  \item{senAUC}{The AUC of the curve given the specificities and the mean sensitivities}
%% ...
}


%%\examples{
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Inspection}
