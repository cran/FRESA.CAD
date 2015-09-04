\name{medianPredict}
\alias{medianPredict}
\title{The median prediction from a list of models}
\description{
	Given a list of model formulas, this function will train such models and return the median prediction on a test data set.
	It also provides a \emph{k}-nearest neighbours (KNN) prediction using the features listed in such models.
}
\usage{
	medianPredict(formulaList,
	              trainData,
	              testData = NULL, 
	              predictType = c("prob", "linear"),
	              type = c("LOGIT", "LM", "COX"),
	              Outcome = "CLASS",
	              nk=0,
	              ...)
}
\arguments{
	\item{formulaList}{
		A list made of objects of class \code{formula}, each representing a model formula to be fitted and predicted with
	}
	\item{trainData}{
		A data frame with the data to train the model, where all variables are stored in different columns
	}
	\item{testData}{
		A data frame similar to \code{trainData}, but with the data set to be predicted. If \code{NULL}, \code{trainData} will be used
	}
	\item{predictType}{
		Prediction type: Probability ("prob") or linear predictor ("linear")
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{nk}{
		The number of neighbours used to generate the KNN classification. If zero, \emph{k} is set to the square root of the number of cases. If less than zero, it will not perform the KNN classification
	}
	\item{...}{
		Additional parameters for fitting a \code{glm} object
	}
}
\value{
	\item{medianPredict}{
		A vector with the median prediction for the \code{testData} data set, using the models from \code{formulaList}
	}
	\item{medianKNNPredict}{
		A vector with the median prediction for the \code{testData} data set, using the KNN models
	}
	\item{predictions}{
		A matrix, where each column represents the predictions made with each model from \code{formulaList}
	}
	\item{KNNpredictions}{
		A matrix, where each column represents the predictions made with a different KNN model
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\examples{
	\dontrun{
	# Start the graphics device driver to save all plots in a pdf format
	pdf(file = "Example.pdf")
	# Get the stage C prostate cancer data from the rpart package
	library(rpart)
	data(stagec)
	# Split the stages into several columns
	dataCancer <- cbind(stagec[,c(1:3,5:6)],
	                    gleason4 = 1*(stagec[,7] == 4),
	                    gleason5 = 1*(stagec[,7] == 5),
	                    gleason6 = 1*(stagec[,7] == 6),
	                    gleason7 = 1*(stagec[,7] == 7),
	                    gleason8 = 1*(stagec[,7] == 8),
	                    gleason910 = 1*(stagec[,7] >= 9),
	                    eet = 1*(stagec[,4] == 2),
	                    diploid = 1*(stagec[,8] == "diploid"),
	                    tetraploid = 1*(stagec[,8] == "tetraploid"),
	                    notAneuploid = 1-1*(stagec[,8] == "aneuploid"))
	# Remove the incomplete cases
	dataCancer <- dataCancer[complete.cases(dataCancer),]
	# Load a pre-stablished data frame with the names and descriptions of all variables
	data(cancerVarNames)
	# Rank the variables:
	# - Analyzing the raw data
	# - According to the zIDI
	rankedDataCancer <- univariateRankVariables(variableList = cancerVarNames,
	                                            formula = "Surv(pgtime, pgstat) ~ 1",
	                                            Outcome = "pgstat",
	                                            data = dataCancer,
	                                            categorizationType = "Raw",
	                                            type = "COX",
	                                            rankingTest = "zIDI",
	                                            description = "Description")
	# Get a Cox proportional hazards model using:
	# - The top 7 ranked variables
	# - 10 bootstrap loops in the feature selection procedure
	# - The zIDI as the feature inclusion criterion
	# - 5 bootstrap loops in the backward elimination procedure
	# - A 5-fold cross-validation in the feature selection, 
	#            update, and backward elimination procedures
	# - A 10-fold cross-validation in the model validation procedure
	# - First order interactions in the update procedure
	cancerModel <- crossValidationFeatureSelection_Bin(size = 7,
	                                               loops = 10,
	                                               Outcome = "pgstat",
	                                               timeOutcome = "pgtime",
	                                               variableList = rankedDataCancer,
	                                               data = dataCancer,
	                                               type = "COX",
	                                               selectionType = "zIDI",
	                                               elimination.bootstrap.steps = 5,
	                                               trainRepetition = 5,
	                                               CVfolds = 10,
	                                               interaction = c(1,2))
	# Get the median prediction:
	# - Without an independent test set
	# - Without a KNN classification
	mp <- medianPredict(formulaList = cancerModel$formula.list,
	                    trainData = dataCancer,
	                    predictType = "prob",
	                    type = "COX",
	                    Outcome = "pgstat",
	                    nk=0)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
