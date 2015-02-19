\name{getKNNpredictionFromFormula}
\alias{getKNNpredictionFromFormula}
\title{Predict classification using KNN}
\description{
	This function will return the classification of the samples of a test set using a \emph{k}-nearest neighbors (KNN) algorithm with euclidean distances, given a formula and a train set.
}
\usage{
	getKNNpredictionFromFormula(model.formula,
	                            trainData,
	                            testData,
	                            Outcome = "CLASS",
	                            nk = 3)
}
\arguments{
	\item{model.formula}{
		An object of class \code{formula} with the formula to be used
	}
	\item{trainData}{
		A data frame with the data to train the model, where all variables are stored in different columns
	}
	\item{testData}{
		A data frame similar to \code{trainData}, but with the data set to be predicted
	}
	\item{Outcome}{
		The name of the column in \code{trainData} that stores the variable to be predicted by the model
	}
	\item{nk}{
		The number of neighbors used to generate the KNN classification
	}
}
\value{
	\item{prediction}{
		A vector with the predicted outcome for the \code{testData} data set
	}
	\item{prob}{
		The proportion of \emph{k} neighbours that predicted the class to be the one being reported in \code{prediction}
	}
	\item{binProb}{
		The proportion of \emph{k} neighbours that predicted the class of the outcome to be equal to 1
	}
	\item{featureList}{
		A vector with the names of the features used by the KNN procedure
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{predictForFresa},
				 \link{knn}}}
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
	# - Age as a covariate
	# - zIDI as the feature inclusion criterion
	cancerModel <- ReclassificationFRESA.Model(loops = 10,
	                                           covariates = "1 + age",
	                                           Outcome = "pgstat",
	                                           variableList = cancerVarNames,
	                                           data = trainDataCancer,
	                                           type = "COX",
	                                           timeOutcome = "pgtime",
	                                           selectionType = "zIDI")
	# Predict the outcome of the test data sample using KNN
	KNNPrediction <- getKNNpredictionFromFormula(model.formula = cancerModel$formula,
	                                             trainData = trainDataCancer,
	                                             testData = testDataCancer,
	                                             Outcome = "pgstat",
	                                             nk = 5)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
