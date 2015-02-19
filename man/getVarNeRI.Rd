\name{getVarNeRI}
\alias{getVarNeRI}
\title{Analysis of the effect of each term of a linear regression model by analyzing its residuals}
\description{
	This function provides an analysis of the effect of each model term by comparing the residuals of the full model and the model without each term.
	The model is fitted using the train data set, but analysis of residual improvement is done on the train and test data sets.
	Residuals are compared by a paired \emph{t}-test, a paired Wilcoxon rank-sum test, a binomial sign test and the \emph{F}-test on residual variance.
	Additionally, the net residual improvement (NeRI) of each model term is reported.
}
\usage{
	getVarNeRI(object,
	           data,
	           Outcome = "Class",
	           type = c("LM", "LOGIT", "COX"),
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
	\item{tP.value}{
		A vector in which each element represents the single sided \emph{p}-value of the paired \emph{t}-test comparing the absolute values of the residuals obtained with the full model and the model without one term
	}
	\item{BinP.value}{
		A vector in which each element represents the \emph{p}-value associated with a significant improvement in residuals according to the binomial sign test
	}
	\item{WilcoxP.value}{
		A vector in which each element represents the single sided \emph{p}-value of the Wilcoxon rank-sum test comparing the absolute values of the residuals obtained with the full model and the model without one term
	}
	\item{FP.value}{
		A vector in which each element represents the single sided \emph{p}-value of the \emph{F}-test comparing the residual variances of the residuals obtained with the full model and the model without one term
	}
	\item{NeRIs}{
		A vector in which each element represents the net residual improvement between the full model and the model without one term
	}
	\item{testData.tP.value}{
		A vector similar to \code{tP.value}, where values were estimated in \code{testdata}
	}
	\item{testData.BinP.value}{
		A vector similar to \code{BinP.value}, where values were estimated in \code{testdata}
	}
	\item{testData.WilcoxP.value}{
		A vector similar to \code{WilcoxP.value}, where values were estimated in \code{testdata}
	}
	\item{testData.FP.value}{
		A vector similar to \code{FP.value}, where values were estimated in \code{testdata}
	}
	\item{testData.NeRIs}{
		A vector similar to \code{NeRIs}, where values were estimated in \code{testdata}
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{getVarReclassification}}}
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
	# - The Wilcoxon rank-sum test as the feature inclusion criterion
	cancerModel <- NeRIBasedFRESA.Model(loops = 10,
	                                    covariates = "1 + age",
	                                    Outcome = "pgstat",
	                                    variableList = cancerVarNames,
	                                    data = trainDataCancer,
	                                    type = "COX",
	                                    testType= "Wilcox",
	                                    timeOutcome = "pgtime")
	# Get the NeRI of each model term in the train data set and in the independent data set
	cancerModelNeRI <- getVarNeRI(object = cancerModel$final.model,
	                              data = testDataCancer,
	                              Outcome = "pgstat",
	                              type = "COX")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
