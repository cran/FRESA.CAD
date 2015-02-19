\name{summaryReport}
\alias{summaryReport}
\title{Report the univariate analysis, the cross-validation analysis and the correlation analysis}
\description{
	This function takes the variables of the cross-validation analysis and extracts the results from the univariate and correlation analyses.
	Then, it prints the cross-validation results, the univariate analysis results, and the correlated variables.
	As output, it returns a list of each one of these results.
}
\usage{
	summaryReport(univariateObject,
	              summaryBootstrap,
	              listOfCorrelatedVariables = NULL,
	              digits = 2)
}
\arguments{
	\item{univariateObject}{
		A data frame that contains the results of the \code{univariateRankVariables} function
	}
	\item{summaryBootstrap}{
		A list that contains the results of the \code{summary.bootstrapValidation} function
	}
	\item{listOfCorrelatedVariables}{
		A matrix that contains the \code{correlated.variables} value from the results obtained with the \code{listTopCorrelatedVariables} function
	}
	\item{digits}{
		The number of significant digits to be used in the print function
	}
}
\value{
	\item{performance.table}{
		A matrix with the tabulated results of the blind test accuracy, sensitivity, specificities, and area under the ROC curve
	}
	\item{coefStats}{
		A data frame that lists all the model features along with its univariate statistics and bootstrapped coefficients
	}
	\item{cor.varibles}{
		A matrix that lists all the features that are correlated to the model variables
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{summary.bootstrapValidation}}}
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
	# Perform a univariate analysis
	rankedDataCancer <- univariateRankVariables(variableList = cancerVarNames,
	                                           formula = "Surv(pgtime, pgstat) ~ 1",
	                                           Outcome = "pgstat",
	                                           data = dataCancer, 
	                                           categorizationType = "Raw", 
	                                           type = "COX", 
	                                           rankingTest = "zIDI",
	                                           description = "Description")
	# Get the variables that have a correlation coefficient 
	# larger than 0.65 at a p-value of 0.05
	cor <- listTopCorrelatedVariables(variableList = cancerVarNames,
	                                  data = dataCancer,
	                                  pvalue = 0.05,
	                                  corthreshold = 0.65,
	                                  method = "pearson")
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
	# Get the summary report
	sumReport <- summaryReport(univariateObject = rankedDataCancer,
	                           summaryBootstrap = sumBootCancerModel,
	                           listOfCorrelatedVariables = cor$correlated.variables)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
