\name{reportEquivalentVariables}
\alias{reportEquivalentVariables}
\title{Report the set of variables that will perform an equivalent IDI discriminant function}
\description{
	Given a model, this function will report a data frame with all the variables that may be interchanged in the model without affecting its classification performance.
	For each variable in the model, this function will loop all candidate variables and report all of which result in an equivalent or better zIDI than the original model.
}
\usage{
	reportEquivalentVariables(object,
	                          pvalue = 0.05,
	                          data,
	                          variableList,
	                          Outcome = "Class", 
	                          type = c("LOGIT", "LM", "COX"),
	                          eqFrac = 0.9,
	                          description = ".")
}
\arguments{
	\item{object}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the IDI , allowed for a pair of variables to be considered equivalent
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{eqFrac}{
		A fraction to which the \emph{z}-score will be relaxed, for a pair of variables to be considered equivalent
	}
	\item{description}{
		The name of the column in \code{variableList} that stores the variable description
	}
}
\value{
	A data frame with three columns. The first column is the original variable of the model. 
	The second column lists all variables that, if interchanged, will not statistically affect the performance of the model.
	The third column lists the corresponding \emph{z}-scores of the IDI for each equivalent variable.
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
	# Get a Cox proportional hazards model using:
	# - 10 bootstrap loops
	# - zIDI as the feature inclusion criterion
	cancerModel <- ForwardSelection.Model.Bin(loops = 10,
	                                           Outcome = "pgstat",
	                                           variableList = cancerVarNames,
	                                           data = dataCancer,
	                                           type = "COX",
	                                           timeOutcome = "pgtime",
	                                           selectionType = "zIDI")
	# Get a data frame with variables that could be interchanged:
	# - Relaxing by a factor of 0.7 the z-score of the IDI
	eqVars <- reportEquivalentVariables(object = cancerModel$final.model,
	                                    data = dataCancer,
	                                    variableList = cancerVarNames,
	                                    Outcome = "pgstat", 
	                                    type = "COX",
	                                    eqFrac = 0.7,
	                                    description = "Description")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
