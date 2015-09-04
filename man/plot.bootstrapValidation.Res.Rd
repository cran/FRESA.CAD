\name{plot.bootstrapValidation_Res}
\alias{plot.bootstrapValidation_Res}
\title{Plot ROC curves of bootstrap results}
\description{
	This function plots ROC curves and a Kaplan-Meier curve (when fitting a Cox proportional hazards regression model) of a bootstrapped model.
}
\usage{
	\method{plot}{bootstrapValidation_Res}(x,
	     xlab = "Years",
	     ylab = "Survival",
	     ...)
}
\arguments{
	\item{x}{
		A \code{bootstrapValidation_Res} object
	}
	\item{xlab}{
		The label of the \emph{x}-axis
	}
	\item{ylab}{
		The label of the \emph{y}-axis
	}
	\item{...}{
		Additional parameters for the plot
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{plot.bootstrapValidation_Bin}}}
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
	# - The Wilcoxon rank-sum test as the feature inclusion criterion
	cancerModel <- ForwardSelection.Model.Res(loops = 10,
	                                    Outcome = "pgstat",
	                                    variableList = cancerVarNames,
	                                    data = dataCancer,
	                                    type = "COX",
	                                    testType= "Wilcox",
	                                    timeOutcome = "pgtime")
	# Bootstrap the parameters of the previous model
	bootCancerModel <- bootstrapValidation_Res(loops = 50,
	                                           model.formula = cancerModel$formula,
	                                           Outcome = "pgstat",
	                                           data = dataCancer,
	                                           type = "COX")
	# Plot the bootstrap results
	plot(x = bootCancerModel)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
