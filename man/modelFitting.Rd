\name{modelFitting}
\alias{modelFitting}
\title{Fit a model to the data}
\description{
	This function fits a linear, logistic, or Cox proportional hazards regression model to given data
}
\usage{
	modelFitting(model.formula,
	             data,
	             type = c("LOGIT", "LM", "COX"),
	             fast=FALSE,
	              ...)
}
\arguments{
	\item{model.formula}{
		An object of class \code{formula} with the formula to be used
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{fast}{
		if true it will perform a fast fitting.
	}
	\item{...}{
		Additional parameters for fitting a default \code{glm} object
	}
}
\value{
	A fitted model of the type defined in \code{type}
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
	# Create a formula of a Cox proportional hazards model using all variables
	allVars <- formula("Surv(pgtime, pgstat)  ~ 1 +
	                                            age +
	                                            g2 +
	                                            grade +
	                                            gleason4 +
	                                            gleason5 +
	                                            gleason6 +
	                                            gleason7 +
	                                            gleason8 +
	                                            gleason910 +
	                                            eet +
	                                            diploid +
	                                            tetraploid +
	                                            notAneuploid")
	# Fit the model to the dataCancer
	allVarsFit <- modelFitting(model.formula = allVars,
	                           data = dataCancer,
	                           type = "COX")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
