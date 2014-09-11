\name{NeRIBasedFRESA.Model}
\alias{NeRIBasedFRESA.Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Net-Residual-Improvement feature selection procedure for regression models (Linear, Logistic or COX) from a set of candidate variables.  
%%  ~~function to do ... ~~
}
\description{
This function preform a Bootstrap sampling to rank the most frequent variables that 
statistically aid the regression models to move towards a better residual. 
After the frequency rank, the function uses a forward selection procedure to create a final linear model 
with terms that provide a significant improvement in residuals.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NeRIBasedFRESA.Model(size = 100, 
fraction = 1, 
pvalue = 0.05, 
loops = 100, 
covariates = "1", 
Outcome, 
variableList, 
dataframe, 
maxTrainModelSize = 10, 
type = c("LM", "LOGIT", "COX"), 
testType=c("Binomial","Wilcox","tStudent","Ftest"),
timeOutcome = "Time",
loop.threshold=20,
interaction = 1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
  The number of candidate variables to be tested.

}
  \item{fraction}{
  The fraction of samples (with replacement) to be used for training and blind test.

}
  \item{pvalue}{
  The critical p.value to control the false selection rate (FSR) 
}
  \item{loops}{
  The number of bootstraps steps to be performed. 

}
  \item{covariates}{
  A string of variables that will be included in the models

}
  \item{Outcome}{
  The name of the outcome variable
%%     ~~Describe \code{Outcome} here~~
}
  \item{variableList}{
  A list with the name of the candidate variables to be tested. 

}
  \item{dataframe}{
  The training data frame whose column include the outcome, the covariates and the variable list. The dataframe rows have all the subject samples

}
  \item{maxTrainModelSize}{
  The maximum size of the model

}
  \item{type}{
  The type of fitting (LOGIT,LM, COX): LOGIT: Logistic, LM Linear, Or COX

}
  \item{testType}{
  The pvalue that will be returned by the improved residuals function
}
  \item{timeOutcome}{
  The name of the time to event variable required for COX fitting.

}
  \item{loop.threshold}{
  The loop.threshold indicates how many bootstrapping cycles will be performed before using the
  selection frequency to filter for candidate features to be included in the model. 
  If the selection frequency of a feature is less than 1/(2*loop.threshold), then it will be removed from the variable search list
}
  \item{interaction}{
  the order of the linear model. Default = 1. Max =2. If set to 2, it will automatically search for pair-wise feature interactions.
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{
%%  If it is a LIST, use
  \item{final.model}{The final model after forward selection}
  \item{var.names}{The list of features that were included in the model}
  \item{formula}{The model formula}
  \item{ranked.var}{The table that contains the feature frequency}
  \item{z.iri}{A vector with the z-value of the residual improvement of each included feature}
  \item{formulas.list}{The list of formulas that were created after each bootstrap cycle}
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

\examples{
## Get the stage C prostate cancer data

require(rpart)
data(stagec)
## Split the stages into several columns
dataCancer <- cbind(stagec[,1:3],stagec[,5:6],1*(stagec[,7] == 4),
1*(stagec[,7] == 5),1*(stagec[,7] == 6),1*(stagec[,7] == 7),
1*(stagec[,7] == 8),1*(stagec[,7] >= 9),1*(stagec[,4] == 2),
1*(stagec[,8]=='diploid'),1*(stagec[,8]=='tetraploid'),1-1*(stagec[,8]=='aneuploid'))
## Assign the column names
colnames(dataCancer) <- c("pgtime","pgstat","age","g2","Grade",
"gleason4","gleason5","gleason6","gleason7","gleason8",
"gleason910","eet","diploid","tetraploid","Notaneuploid")
## Lets remove the incomplete cases
dataCancer <- dataCancer[complete.cases(dataCancer),]
## Load the list of variables to be searched
data(CancerVarNames)
\dontrun{
## Rank the COX model variables according to the IDI
rankedCancer <- univariateRankVariables(CancerVarNames,"Surv(pgtime,pgstat) ~ 1",
 "pgstat",dataCancer,"Raw","COX","zIDI",description="Description");
## Get the bootstrapped model that reduces the  Martingale  residuals of the prognosis
CancerModel <- NeRIBasedFRESA.Model(20,1.0,0.05,1,'1','pgstat',
rankedCancer,dataCancer,type="COX",timeOutcome="pgtime")

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Generation}
