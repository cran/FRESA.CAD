\name{featureAdjustment}
\alias{featureAdjustment}
\title{Adjust each listed variable to the provided set of covariates}
\description{
	This function fits the candidate variables to the provided model,for each strata, on a control population.
	If the variance of the residual (the fitted observation minus the real observation) is reduced significantly, then, such residual is used in the resulting data frame.
	Otherwise, the control mean is subtracted to the observation.
}
\usage{
	featureAdjustment(variableList,
	                  baseModel,
	                  strata = "NA",
	                  data,
	                  referenceframe,
	                  type = c("LM", "GLS"),
	                  pvalue = 0.05,
	                  correlationGroup = "ID")
}
\arguments{
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{baseModel}{
		A string of the type "1 + var1 + var2" that defines the model to which variables will be fitted
	}
	\item{strata}{
		The name of the column in \code{data} that stores the variable that will be used to stratify the model
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{referenceframe}{
		A data frame similar to \code{data}, but with only the control population
	}
	\item{type}{
		Fit type: linear fitting ("LM"), or generalized least squares fitting ("GLS")
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the \emph{F}-test, for the model to be allowed to reduce variability
	}
	\item{correlationGroup}{
		The name of the column in \code{data} that stores the variable to be used to group the data (only needed if \code{type} defined as "GLS")
	}
}
\value{
	A data frame, where each input observation has been adjusted from \code{data} at each \code{strata}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\note{This function prints the residuals and the \emph{F}-statistic for all candidate variables}
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
	# Generate a reference frame
	controls <- dataCancer[which(dataCancer$pgstat == 0),]
	# Load a pre-stablished data frame with the names and descriptions of all variables
	data(cancerVarNames)
	# Adjust the g2 variable to age
	adjDataCancer<-featureAdjustment(variableList = cancerVarNames[2,],
	                                 baseModel = "1 + age",
	                                 data = dataCancer,
	                                 referenceframe = controls,
	                                 type = "LM")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Data_Conditioning}