\name{bootstrapValidation}
\alias{bootstrapValidation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap-validation of Binary Classification models
%%  ~~function to do ... ~~
}
\description{
This function bootstrap the model N times to estimate for each variable the empirical distribution of model coefficients along with their IDI and NRI.
At each bootstrap the non-observed data is predicted by the trained model. The method keeps track of prediction and plots the final Blinded ROC, the blind: Accuracy, Sensitivity and Specificity 
contrasted with the bootstrapped trained distributions. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapValidation(fraction = 1, loops = 200, model.formula, 
Outcome, dataframe, type = c("LM", "LOGIT", "COX"),plots=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fraction}{
  The fraction of data to be sampled with replacement (bootstrap sampling)
}
  \item{loops}{
  The number of bootstrapped steps
}
  \item{model.formula}{
  The formula to be bootstrapped and validated. It has to be either a COX object or a LM model.
%%     ~~Describe \code{model.formula} here~~
}
  \item{Outcome}{
  The outcome variable to which the accuracy, sensitivity and specificity will be computed.
}
  \item{dataframe}{
  The training data frame 
}
  \item{type}{
  The type of fit: LM", "LOGIT", "COX"
}
  \item{plots}{
  Set to FALSE to avoid showing the density plots 
}
}
\details{
The bootstrap validation will estimate the CI of the model coefficients, the NRI and IDI values only if the sample fraction is equal to one.
The non-sampled values will be used to do an estimation of the blind Accuracy, Sensitivity and Specificity.  
This function will also provides plots to monitor the evolution of the bootstrap procedure. The plot will show the bootstrapped histogram of accuracy, sensitivity and specificity. 
The plots also will show the density distribution of the blinded test. 
At the end the function will plot the Train ROC, the Blind ROC, the Final Blind: Accuracy, Sensitivity and Specificity will be 
drawn along the continuous cumulative distribution of the trained results.
%%  ~~ If necessary, more details than the description above ~~
}
\value{

%%  If it is a LIST, use
  \item{bin.Outcome}{A vector of the binary outcome}
  \item{blind.accuracy}{The blind accuracy of the model}
  \item{blind.sensitivity}{The blind sensitivity of the model}
  \item{blind.specificity}{The blind specificity of the model}
  \item{train.ROCACU}{The train ROC AUC}
  \item{blind.ROCAUC}{The blind ROC AUC}
  \item{boot.ROCAUC}{The ROC AUC using the mean of the bootstrapped coefficients}
  \item{fraction}{The sampled fraction}
  \item{loops}{The number of bootstrapped loops}
  \item{base.Accuracy}{The accuracy of the original model}
  \item{base.sensitivity}{The sensitivity of the original model}
  \item{base.specificity}{The specificity of the original model}
  \item{accuracy}{A vector with the bootstrapped test accuracies}
  \item{sensitivities}{A vector with the bootstrapped test sensitivities}
  \item{specificities}{A vector with the bootstrapped test specificities}
  \item{train.accuracy}{A vector with the bootstrapped train accuracies}
  \item{train.sensitivity}{A vector with the bootstrapped train sensitivities}
  \item{train.specificity}{A vector with the bootstrapped train specificities}
  \item{s.coef}{a data matrix with the train coefficients. Each row contains a bootstrapped sample of the coefficients}
  \item{boot.model}{The model whose coefficients are the mean of the bootstrapped models}
  \item{mboot.model}{The model whose coefficients are the IDI weighted mean of the bootstrapped models}
  \item{boot.accuracy}{The estimated test accuracy}
  \item{boot.sensitivity}{The estimated test sensitivity}
  \item{boot.specificity}{The estimated test specificity}
  \item{z.NRIs}{A data matrix of the estimated zNRI for each model coefficient. Each row is a bootstrap sample}
  \item{z.IDIs}{A data matrix of the estimated zIDI for each model coefficient. Each row is a bootstrap sample }
  \item{test.z.NRIs}{A data matrix of the estimated test zNRI for each model coefficient. Each row is a bootstrap sample}
  \item{test.z.IDIs}{A data matrix of the estimated test zIDI for each model coefficient. Each row is a bootstrap sample }
  \item{NRIs}{A data matrix of the estimated NRI for each model coefficient. Each row is a bootstrap sample}
  \item{IDIs}{A data matrix of the estimated IDI for each model coefficient. Each row is a bootstrap sample }
  \item{blind.test.outcome}{A vector that contains all the individual outcomes used to validate the model}
  \item{blind.test.prediction}{A vector that contains all the individual predictions used to validate the model}
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
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Diagnosis}
