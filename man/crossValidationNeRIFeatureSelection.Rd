\name{crossValidationNeRIFeatureSelection}


\alias{crossValidationNeRIFeatureSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IDI-based selection of a linear classification (Logistic or COX) model from a set of candidate variables.  
%%  ~~function to do ... ~~
}
\description{
This function preform the cross-validation analysis of the NeRI based feature selection algorithm.
Internally is composed of a NeRI selection, followed by an update procedure and finalized by a backwards elimination procedure.
The number of cross-validation folds in given by the size of the training set. The user can control how many steps the training and testing procedure will be repeated.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
crossValidationNeRIFeatureSelection(
size=10,
fraction=1.0,
pvalue=0.05,
loops=100,
covariates="1",
Outcome,
variableList,
dataframe,
maxTrainModelSize=10,
type=c("LM","LOGIT","COX"),
testType=c("Binomial","Wilcox","tStudent","Ftest"),
timeOutcome="Time",
loop.threshold=10,
startOffset=0,
trainFraction=0.67,
trainRepetition=9,
elimination.pValue=0.05,
setIntersect=1,
interaction=c(1,1),
update.pvalue=c(0.05,0.05),
backBootLoops=25
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
  The number of candidate variables to be tested.

}
  \item{fraction}{
  The fraction of samples used at the training and the independent testing .

}
  \item{pvalue}{
  The threshold used to insert a variable in the model. the p value of the N(z.IDI,0,1) has to be lower that pvalue in order to be included in the model.

}
  \item{loops}{
  The number of bootstraps 

}
  \item{covariates}{
  A string of variables that always will be included in the models

}
  \item{Outcome}{
  The variable that indicates the outcome

}
  \item{variableList}{
  The candidate variables to be tested. They usualy are the output of a univariate ranking procedure

}
  \item{dataframe}{
  The data frame with all variables and co-variables in columns and row samples

}
  \item{maxTrainModelSize}{
  The maximum size of models to be tested

}
  \item{type}{
  The type of fitting (LOGIT,LM, COX): LOGIT: Logistic, LM: Lineal, Or COX

}
  \item{testType}{
  The type of test for the improvedResiduals function

}
  \item{timeOutcome}{
  The variable that has the time to event for COX fitting.

}
  \item{loop.threshold}{
  it defines the number of loops to be used to check for the frequency of variables that will continue the bootstrapping procedure
	}
  \item{startOffset}{
  the number of items to be keep in the backwards elimination procedure
  }
  \item{trainFraction}{
  the fraction of subjects that will be used in each training step of the cross validation fold
  }
  \item{setIntersect}{
  Set the intersect for the update model
  }
  \item{trainRepetition}{
  the number of cross validation steps. For a complete cross validation it should be equal to 1/crossFraction
  }
  \item{elimination.pValue}{
  the maximum p.value allowed for the z.IDI of a model variable.
  }
  \item{interaction}{
  the order of the linear model. Default = 1. Max =2. If set to 2, it will search for pair-wise feature interactions.
	}
  \item{update.pvalue}{
  The p.value for model update
	}
  \item{backBootLoops}{
  The number of bootstrap steps to be performed
	}

}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{

%%  If it is a LIST, use
  \item{formula.list}{A list with the formulas found at each CV step}
  \item{Models.testPrediction}{A matrix with the test results for each CV fold}
  \item{FullModel.testPrediction}{A matrix with the test results for the model that observed all the data}
  \item{FullBackModel}{The final model that observed all the data}
  \item{FullBaggingModel}{The full model after bagging}
  \item{FullUpdatedModel}{The model after the update procedure}
  \item{blindRMS}{The RMS of the blind folds}
  \item{blindPearson}{The blind Pearson correlation}
  \item{FullblindRMS}{The re-test RMS value of the full model}
  \item{FullblindPearson}{The re-test Pearson value of the full model}
  \item{FullblindSpearman}{The re-test Spearman value of the full model}
  \item{TrainRMS}{The train RMS values}
  \item{TrainSpearman}{The train Speaman values}
  \item{TrainPearson}{The train Pearson}
  \item{TrainRMS.Full}{The re-train  RMS values}
  \item{TrainSpearman.FULL}{The re-train spearman}
  \item{TrainPearson.FULL}{The re-train pearson}
  \item{BlindRMSAtFold}{Per fold blind per fold RMS values}
  \item{BlindRMSAtFold.FULL}{Per fold re-test RMS values}
	\item{fullenet}{The elastic-net cross-validation object}
	\item{enet.testPredictions}{The elastic-net fold predictions}
	\item{enetVariables}{The elastic-net variables}
  
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
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Generation}
