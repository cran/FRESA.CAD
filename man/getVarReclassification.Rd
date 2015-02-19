\name{getVarReclassification}
\alias{getVarReclassification}
\title{Analysis of the effect of each term of a binary classification model by analyzing its reclassification performance}
\description{
	This function provides an analysis of the effect of each model term by comparing the binary classification performance between the full model and the model without each term.
	The model is fitted using the train data set, but probabilities are predicted for the train and test data sets.
	Reclassification improvement is evaluated using the \code{improveProb} function (\code{Hmisc} package).
	Additionally, the integrated discrimination improvement (IDI) and the net reclassification improvement (NRI) of each model term are reported.
}
\usage{
	getVarReclassification(object,
	                       data,
	                       Outcome = "Class", 
	                       type = c("LOGIT", "LM", "COX"),
	                       testData = NULL)
}
\arguments{
	\item{object}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{testData}{
		A data frame similar to \code{data}, but with a data set to be independently tested. If \code{NULL}, \code{data} will be used.
	}
}
\value{
	\item{z.IDIs}{
		A vector in which each term represents the \emph{z}-score of the IDI obtained with the full model and the model without one term
	}
	\item{z.NRIs}{
		A vector in which each term represents the \emph{z}-score of the NRI obtained with the full model and the model without one term
	}
	\item{IDIs}{
		A vector in which each term represents the IDI obtained with the full model and the model without one term
	}
	\item{NRIs}{
		A vector in which each term represents the NRI obtained with the full model and the model without one term
	}
	\item{testData.z.IDIs}{
		A vector similar to \code{z.IDIs}, where values were estimated in \code{testdata}
	}
	\item{testData.z.NRIs}{
		A vector similar to \code{z.NRIs}, where values were estimated in \code{testdata}
	}
	\item{testData.IDIs}{
		A vector similar to \code{IDIs}, where values were estimated in \code{testdata}
	}
	\item{testData.NRIs}{
		A vector similar to \code{NRIs}, where values were estimated in \code{testdata}
	}
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{getVarNeRI}}}
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
	# Get the IDI and NRI of each model term in the train data 
	# set and in the independent data set
	cancerModelRec <- getVarReclassification(object = cancerModel$final.model,
	                                         data = trainDataCancer,
	                                         Outcome = "pgstat",
	                                         type = "COX",
	                                         testData = testDataCancer)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}