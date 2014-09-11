\name{getVarReclassification}
\alias{getVarReclassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analysis of the effect of each model-term on the binary reclassification performance
%%  ~~function to do ... ~~
}
\description{
Each model term is analysed by comparing the binary-classification performance between the full and reduced model.
The model is fitted on the provided data set and the predicted probability of the training outcome and the independent-test outcome is
evaluated by the improveProb function of the HMISC package.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getVarReclassification(object, dataframe, Outcome = "Class", 
type = c("LOGIT", "LM", "COX"),
independentFrame=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  The model to be tested
}
  \item{dataframe}{
  The data frame to fit

}
  \item{Outcome}{
  The Binary outcome to be tested
}
  \item{type}{
  The type of fit : COX, LOGIT or LM


}
  \item{independentFrame}{
  an independent frame used for the estimation of the NRI/IDI improved probability
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{

%%  If it is a LIST, use
  
  \item{z.IDIs }{The vector of z values for each variable IDI estimated with the independent data}
  \item{z.NRIs }{The vector of z values for each variable NRI estimated with the independent data}
  \item{IDIs }{The vector of IDI for each variable using the independent data}
  \item{NRIs }{The vector of IDI for each variable using the independent data}
  \item{tz.IDIs }{The vector of z values for each variable IDI estimated with the train data}
  \item{tz.NRIs}{The vector of z values for each variable NRI estimated with the train data}
  \item{tIDIs}{The vector of IDI for each variable using the train data}
  \item{tNRIs }{The vector of NRI for each variable using the train data}
%% ...
}
\references{
Pencina, Michael J., Ralph B. D'Agostino, and Ramachandran S. Vasan. 
"Evaluating the added predictive ability of a new marker: from area under 
the ROC curve to reclassification and beyond." 
Statistics in medicine 27.2 (2008): 157-172.
}
\author{
Jose G. Tamez-Pena

}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 See Also as \code{improveProb} of the HMISC package 
}

%%\examples{
%%
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Inspection}
