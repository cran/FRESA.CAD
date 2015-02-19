\name{summary.bootstrapValidation}
\alias{summary}
\alias{summary.bootstrapValidation}
\title{Generate a report of the results obtained using the bootstrapValidation function}
\description{
	This function prints two tables describing the results of the bootstrap-based validation of binary classification models.
	The first table reports the accuracy, sensitivity, specificity and area under the ROC curve (AUC) of the train and test data set, along with their confidence intervals.
	The second table reports the model coefficients and their corresponding integrated discrimination improvement (IDI) and net reclassification improvement (NRI) values.
}
\usage{
	\method{summary}{bootstrapValidation}(object,
	        ...)
}
\arguments{
	\item{object}{
		An object of class \code{bootstrapValidation}
	}
	\item{...}{
		Additional parameters for the generic \code{summary} function
	}
}
\value{
	\item{performance}{
		A vector describing the results of the bootstrapping procedure
	}
	\item{summary}{
		An object of class \code{summary.lm}, \code{summary.glm}, or \code{summary.coxph} containing a summary of the analyzed model
	}
	\item{coef}{
		A matrix with the coefficients, IDI, NRI, and the 95\% confidence intervals obtained via bootstrapping
	}
	\item{performance.table}{
		A matrix with the tabulated results of the blind test accuracy, sensitivity, specificities, and area under the ROC curve
	}  
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{summaryReport}}}
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
	# Get a Cox proportional hazards model using:
	# - 10 bootstrap loops
	# - Age as a covariate
	# - zIDI as the feature inclusion criterion
	cancerModel <- ReclassificationFRESA.Model(loops = 10,
	                                           covariates = "1 + age",
	                                           Outcome = "pgstat",
	                                           variableList = cancerVarNames,
	                                           data = dataCancer,
	                                           type = "COX",
	                                           timeOutcome = "pgtime",
	                                           selectionType = "zIDI")
	# Validate the previous model:
	# - Using 50 bootstrap loops
	bootCancerModel <- bootstrapValidation(loops = 50,
	                                       model.formula = cancerModel$formula,
	                                       Outcome = "pgstat",
	                                       data = dataCancer,
	                                       type = "COX")
	# Get the summary of the bootstrapped model
	sumBootCancerModel <- summary.bootstrapValidation(object = bootCancerModel)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
