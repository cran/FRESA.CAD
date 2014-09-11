\name{crossValidationFeatureSelection}
\alias{crossValidationFeatureSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
IDI-based selection of a linear classification (Logistic or COX) model from a set of candidate variables.  
%%  ~~function to do ... ~~
}
\description{
This function preform the cross-validation analysis of the IDI based feature selection algorithm.
Internally is composed of a IDI selection, followed by an update procedure and finalized by a bootstrapping backwards elimination procedure.
The number of cross-validation folds in given by the size of the training set. The user can control how many steps the training and testing procedure will be repeated.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
crossValidationFeatureSelection(size=10,
fraction=1.0,
pvalue=0.05,
loops=100,
covariates="1",
Outcome,
variableList,
dataframe,
maxTrainModelSize=10,
type=c("LM","LOGIT","COX"),
timeOutcome="Time",
selectionType=c("zIDI","zNRI"),
loop.threshold=10,
startOffset=0,
backBootLoops=25,
trainFraction=0.67,
trainRepetition=9,
elimination.pValue=0.05,
CVfolds=10,
bootEstimations=25,
interaction=c(1,1),
nk = 0
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
  THe number of bootstraps 

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
  The data

}
  \item{maxTrainModelSize}{
  A threshold for the maximum size of models to be tested

}
  \item{type}{
  The type of fitting (LOGIT,LM, COX): LOGIT: Logistic, LM: Lineal, Or COX

}
  \item{timeOutcome}{
  The name of the time to event variable for COX fitting.

}
  \item{selectionType}{
  Selection criteria for feature inclusion: zIDI or zIRI
  }
  \item{loop.threshold}{
  After loop.threshold cycles no more variables will be added to the frequency histogram. After that point only variables with at least one count will be candidates for the feature selection procedure
	}
  \item{startOffset}{
  Only terms whose index in greater than startOffset can be removed by the backwards elimination procedure
  }
  \item{backBootLoops}{
  The number of bootstraps to be used in the backwards elimination procedure
  }
  \item{trainFraction}{
  the fraction of subjects that will be used in each training step of the cross validation fold
  }
  \item{trainRepetition}{
  the number of cross validation steps. For a complete cross validation it should be at least equal to 1/crossFraction
  }
  \item{elimination.pValue}{
  the maximum term IDI p.value allowed in the formula.
  }
  \item{CVfolds}{
  the number of folds for the final CV of generated models.
  }
  \item{bootEstimations}{
  number of bootstrap to be used to estimate each model coefficients 95% CI. 
  }
  \item{interaction}{
  the order of the factors. Default = 1. Max =2. If set to 2, it will search for pair-wise feature interactions.
}
  \item{nk}{
  the number of neighbours for a KNN classification. 
}

}

\details{
This function for the cross-validation (CV) of the feature selection procedure will produce a set data and plots that can be used to inspect the degree of over-fitting or model shrinkage
present on the selected models. For each generated model bootstrapped data, CV data, and, if possible, retrain data. At each cycle the bootstrapped model will generate a training ROC 
and a test ROC of the final model. At the end of the CV of the selection procedure, a set of three plots may be produced. The first plot will contain the ROC for each blind test.
A second possible plot will contain the ROC plots of each model train and tested in the blind fold. The final plot will have the train ROC plot generated with all the data. The 
ROC plot of the bootstrapped blind test, the ROC plot of the CV test, along of a set of plots of the CV of each model. Finally, this last plot will contain the ROC analysis of the 
mean test data, along the coherence ROC plot. These set of plots may be used to get an overall view of the expected model shrinkage.
Along with the plots the function will be providing the overall performance of the system( Accuracy, sensitivity and specificity). The function also produces a report of the expected
performance of a KNN algorithm trained with the selected features of the model. The KNN test prediction are compared to the predictions generated by the linear models.  
}

\value{

%%  If it is a LIST, use
  \item{formula.list}{A vector with the formulas found at each cycle of the cross-validation procedure}
  \item{Models.testPrediction}{A data frame with the results of the blind test for each model: The unbiased expectation of the model when faced with new data}
  \item{FullModel.testPrediction}{A data frame with the result of the cross-validation procedure of the model found using all available data}
  \item{Models.CVtestPredictions}{The data frame of the last CV test predictions. This should model the FullModel CV blind test}
  \item{TestRetrained.blindPredictions}{If enough samples present in the fold set, this data frame will have the data the get a reliable estimation of the model performance retrained on an independent set of data}
  \item{FullModel.bootstrapped}{The bootstrapped object of the model trained using all available information}
  \item{LastTrainedModel.bootstrapped}{The bootstrapped object of the last trained model}
  \item{Test.accuracy}{The expected blind accuracy of the trained model}
  \item{Test.sensitivity}{The expected blind sensitivity of the trained model }
  \item{Test.specificity}{The expected blind specificity of the trained model }
  \item{Train.correlationsToFull}{The degree of correlation between each model and the full model on the training data sets}
  \item{Blind.correlationsToFull}{The degree of correlation between each model and the full model on the test folds}
  \item{FullModelAtFoldAccuracies}{The test accuracy for the full model at each CV fold}
  \item{FullModelAtFoldSpecificties}{The test specificity for the full model at each CV fold}
  \item{FullModelAtFoldSensitivities}{The test sensitivity for the full model at each CV fold}
  \item{AtCVFoldModelBlindAccuracies}{Blind Accuracy of the selection procedure}
  \item{AtCVFoldModelBlindSpecificities}{Blind Specificity of the selection procedure}
  \item{AtCVFoldModelBlindSensitivities}{Blind Sensitivity of the selection procedure}
  \item{CVROCSensitivities}{The mean of the ROC sensitivities at the given specificities \cr
  of: [0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05]}
  \item{Models.CVblindMeanSensitivites}{The mean ROC sensitivities of the blind folds at specificities \cr
  of: [0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05]}
  \item{varFullSelection}{The model returned by the ReclassificationFRESA.Model using all available data}
  \item{Models.testSensitivities }{Each ROC test sensitivities at the given specificities \cr
  of: [0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05]}
  \item{FullKNN.testPrediction}{The predicted outcome by the KNN procedure at each fold using the features of the full model}
  \item{KNN.testPrediction}{The predicted outcome by the KNN procedure at each fold using the fold model}
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
