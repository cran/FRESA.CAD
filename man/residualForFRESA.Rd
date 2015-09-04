\name{residualForFRESA}
\alias{residualForFRESA}
\title{Return residuals from prediction}
\description{
	Given a model and a new data set, this function will return the residuals of the predicted values.
	When dealing with a Cox proportional hazards regression model, the function will return the Martingale residuals.
}
\usage{
	residualForFRESA(object,
	                 testData,
	                 Outcome,
	                 eta = 0.05)
}
\arguments{
	\item{object}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{testData}{
		A data frame  where all variables are stored in different columns, with the data set to be predicted
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{eta}{
		The weight of the contribution of the Martingale residuals, or 1 - the weight of the contribution of the classification residuals (only needed if \code{object} is of class \code{coxph})
	}
}
\value{
	A vector with the residuals (i.e. the differences between the predicted and the real outcome)
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
	# - The ranked variables
	# - The Wilcoxon rank-sum test as the feature inclusion criterion
	cancerModel <- ForwardSelection.Model.Res(loops = 10,
	                                    Outcome = "pgstat",
	                                    variableList = cancerVarNames,
	                                    data = trainDataCancer,
	                                    type = "COX",
	                                    testType= "Wilcox",
	                                    timeOutcome = "pgtime")
	# Get the residuals of the model
	# - In the test data
	# - Giving the same weight to the Martingale and classification residuals
	cancerModelRes <- residualForFRESA(object = cancerModel$final.model,
	                                   testData = testDataCancer,
	                                   Outcome = "pgstat",
	                                   eta = 0.5)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
