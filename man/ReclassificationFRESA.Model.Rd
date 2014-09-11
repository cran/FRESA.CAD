\name{ReclassificationFRESA.Model}
\alias{ReclassificationFRESA.Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IDI-based selection of a linear classification (Logistic or COX) model from a set of candidate variables.  
%%  ~~function to do ... ~~
}
\description{
This function preforms a bootstrap sampling to rank the variables that statistically improve predictor models using the train/test bootstrapped samples
At each bootstrap step each variable IDI/NRI is computed and the variables with the largest statically significance IDI/NRI is added to the model.
The procedure is repeated at each bootstrap until no more variables can be inserted.
The variables that enter the model are counted and the next bootstrap sample is repeated for N times. 
The frequency of variable-inclusion in the model is returned as well as a model that uses the frequency of inclusion.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ReclassificationFRESA.Model(size = 100, fraction = 1.0, pvalue = 0.05, 
loops = 100, covariates = "1", Outcome, variableList, dataframe, 
maxTrainModelSize = 10, type = c("LM", "LOGIT", "COX"), timeOutcome = "Time",
selectionType=c("zIDI","zNRI"),
loop.threshold=20,interaction=1
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
  The number of candidate variables to be tested.

}
  \item{fraction}{
  The fraction of samples used at the training and the independent testing phase.

}
  \item{pvalue}{
  The threshold used to insert a variable in the model. the p value of the N(z.IDI,0,1) has to be lower that pvalue in order to be included in the model.

}
  \item{loops}{
  The number of bootstraps steps

}
  \item{covariates}{
  A string of variables that always will be included in the models

}
  \item{Outcome}{
  The outcome variable

}
  \item{variableList}{
  The candidate variables to be tested. They usually are the output of a univariate ranking procedure

}
  \item{dataframe}{
  The train data

}
  \item{maxTrainModelSize}{
  A largest size of models to be created

}
  \item{type}{
  The type of fitting (LOGIT,LM, COX): LOGIT: Logistic, LM: Lineal, Or COX

}
  \item{timeOutcome}{
  Time to event variable for COX fitting.

}
  \item{selectionType}{
  Selection criteria for feature inclusion: zIDI or zIRI
}
  \item{loop.threshold}{
  It defines the number of bootstraps loops to be used before start filtering by variable frequency. 
  Once the bootstrap iteration number is greater than loop.threshold, only variables with a bagging frequency greater than 1/(2*loop.threshold) 
  will be explored as model candidates. 
}
  \item{interaction}{
  the order of the factors. Default = 1. Max =2. If set to 2, it will search for first order feature interactions.
}

}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{

%%  If it is a LIST, use
  \item{final.model}{Description of 'comp1'}
  \item{var.names}{Description of 'comp2'}
  \item{formula}{Description of 'comp2'}
  \item{ranked.var}{Description of 'comp2'}
  \item{z.idi}{Description of 'comp2'}
  \item{formula.list}{Description of 'comp2'}
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

%%seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

%%\examples{
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Generation}
