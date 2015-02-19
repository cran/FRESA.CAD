\name{updateModel}
\alias{updateModel}
\title{Update the IDI/NRI-based model using new data or new threshold values}
\description{
	This function will take the frequency-ranked set of variables and will generate a new model with terms that meet either the integrated discrimination improvement (IDI), or the net reclassification improvement (NRI), threshold criteria.
}
\usage{
	updateModel(Outcome,
	            covariates = "1",
	            pvalue = c(0.05, 0.02),
	            VarFrequencyTable, 
	            variableList,
	            data,
	            type = c("LM", "LOGIT", "COX"), 
	            lastTopVariable = 0,
	            timeOutcome = "Time",
	            selectionType = c("zIDI","zNRI"),
	            numberOfModels = 3,
	            interaction = 1,
	            maxTrainModelSize = 0)
}
\arguments{
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{covariates}{
		A string of the type "1 + var1 + var2" that defines which variables will always be included in the models (as covariates)
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to either IDI or NRI, allowed for a term in the model
	}
	\item{VarFrequencyTable}{
		An array with the ranked frequencies of the features, (e.g. the \code{ranked.var} value returned by the \code{ReclassificationFRESA.Model} function)
	}
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{lastTopVariable}{
		The maximum number of variables to be tested
	}
	\item{timeOutcome}{
		The name of the column in \code{data} that stores the time to event (needed only for a Cox proportional hazards regression model fitting)
	}
	\item{selectionType}{
		The type of index to be evaluated by the \code{improveProb} function (\code{Hmisc} package): \emph{z}-score of IDI or of NRI
	}
	\item{numberOfModels}{
		The number of models to be extracted based on the ranked variables
	}
	\item{interaction}{
		Set to either 1 for first order models, or to 2 for second order models
	}
	\item{maxTrainModelSize}{
		Maximum number of terms that can be included in the model
	}
}
\value{
	\item{final.model}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the final model
	}
	\item{var.names}{
		A vector with the names of the features that were included in the final model
	}
	\item{formula}{
		An object of class \code{formula} with the formula used to fit the final model
	}
	\item{z.selectionType}{
		A vector in which each term represents the \emph{z}-score of the index defined in \code{selectionType} obtained with the full model and the model without one term
	}
	\item{loops}{
		The number of loops it took for the model to stabilize
	}
	\item{formula.list}{
		A list containing objects of class \code{formula} with the formulas used to fit the models found at each cycle
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{updateNeRIModel}}}
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
	cancerModel <- ReclassificationFRESA.Model(loops = 10,
	                                           Outcome = "pgstat",
	                                           variableList = cancerVarNames,
	                                           data = dataCancer,
	                                           type = "COX",
	                                           timeOutcome = "pgtime",
	                                           selectionType = "zIDI")
	# Update the model, adding first order interactions
	uCancerModel <- updateModel(Outcome = "pgstat",
	                            VarFrequencyTable = cancerModel$ranked.var,
	                            variableList = cancerVarNames,
	                            data = dataCancer,
	                            type = "COX",
	                            timeOutcome = "pgtime",
	                            interaction = 2)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Generation}