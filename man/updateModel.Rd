\name{updateModel}
\alias{updateModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Updates the IDI/NRI-based generated model using new data or new threshold values.
%%  ~~function to do ... ~~
}
\description{
The update function will take the frequency-ranked set of variables and will generate a new model whose terms meet the IDI/NRI threshold criteria.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
updateModel(Outcome, covariates = "1", pvalue = c(0.05, 0.02), VarFrequencyTable, 
variableListNames, dataframe, type = c("LM", "LOGIT", "COX"), 
lastTopVariable = 0, timeOutcome = "Time",
selectionType=c("zIDI","zNRI"),numberOfModels=3,interaction=1,maxTrainModelSize=0
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Outcome}{
  The binary outcome
}
  \item{covariates}{
  A string of the co-variates that will always be present in the models
}
  \item{pvalue}{
  The p.value to decide if the frequency ranked variable will be included in the model
}
  \item{VarFrequencyTable}{
  A data frame where with an ID of ranked set of variables 
%%     ~~Describe \code{VarFrequencyTable} here~~
}
  \item{variableListNames}{
  The names of the ranked variables
%%     ~~Describe \code{variableListNames} here~~
}
  \item{dataframe}{
  The data frame to be used
}
  \item{type}{
  The fitting type: "LM", "LOGIT" or "COX"
  Specify if the data will be fitted via LM, LOGIT: Logistic or COX models
%%     ~~Describe  \code{type} here~~
}
  \item{lastTopVariable}{
  The maximum number of variables to be tested
%%     ~~Describe \code{lastTopVariable} here~~
}
  \item{timeOutcome}{
  COX fitting only: The variable that contains the time to event information.
}
  \item{selectionType}{
  Selection criteria for feature inclusion: zIDI or zIRI
}
  \item{numberOfModels}{
  set the number of models to be extracted based on the ranked variables
}
  \item{interaction}{
  the order of the linear model. Default = 1. Max =2. If set to 2, it will automatically search for pair-wise feature interactions.
}
  \item{maxTrainModelSize}{
  the maximum number of independent variables in the models
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{

%%  If it is a LIST, use
  \item{final.model}{The Updated Model}
  \item{var.names}{The included features in the model}
  \item{formula}{The model formula}
  \item{z.idi}{the IDI z value for each model term}
  \item{loops}{The number of loops required to stabilize the model}
  \item{formulas}{The list of model formulas}
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
