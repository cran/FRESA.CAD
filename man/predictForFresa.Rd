\name{predictForFresa}
\alias{predictForFresa}
\title{Linear or probabilistic prediction}
\description{
	This function returns the predicted outcome of a specific model.
	The model is used to generate linear predictions.
	The probabilistic values are generated using the logistic transformation on the linear predictors.
}
\usage{
	predictForFresa(object,
	                testData,
	                predictType = c("prob", "linear"))
}
\arguments{
	\item{object}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{testData}{
		A data frame  where all variables are stored in different columns, with the data set to be predicted
	}
	\item{predictType}{
		Prediction type: Probability ("prob") or linear predictor ("linear")
	}
}
\value{
	A vector with the predicted values
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
	# Split the data set into train and test samples
	trainDataCancer <- dataCancer[1:(nrow(dataCancer)/2),]
	testDataCancer <- dataCancer[(nrow(dataCancer)/2+1):nrow(dataCancer),]
	# Get a Cox proportional hazards model using:
	# - 10 bootstrap loops
	# - Train data
	# - zIDI as the feature inclusion criterion
	# - First order interactions
	cancerModel <- ForwardSelection.Model.Bin(loops = 10,
	                                           Outcome = "pgstat",
	                                           variableList = cancerVarNames,
	                                           data = trainDataCancer,
	                                           type = "COX",
	                                           timeOutcome = "pgtime",
	                                           selectionType = "zIDI",
	                                           interaction = 2)
	# Predict the outcome of the test data sample
	predTest <- predictForFresa(object = cancerModel$final.model,
	                            testData = testDataCancer,
	                            predictType = "prob")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
