\name{updateNeRIModel}
\alias{updateNeRIModel}
\title{Update the NeRI-based model using new data or new threshold values}
\description{
	This function will take the frequency-ranked set of variables and will generate a new model with terms that meet the net residual improvement (NeRI) threshold criteria.
}
\usage{
	updateNeRIModel(Outcome, 
	                covariates = "1", 
	                pvalue = c(0.05, 0.02),
	                VarFrequencyTable, 
	                variableList, 
	                data, 
	                type = c("LM", "LOGIT", "COX"),
	                testType=c("Binomial", "Wilcox", "tStudent"), 
	                lastTopVariable = 0, 
	                timeOutcome = "Time",
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
		The maximum \emph{p}-value, associated to the NeRI, allowed for a term in the model
	}
	\item{VarFrequencyTable}{
		An array with the ranked frequencies of the features, (e.g. the \code{ranked.var} value returned by the \code{NeRIBasedFRESA.Model} function)
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
	\item{testType}{
		Type of non-parametric test to be evaluated by the \code{improvedResiduals} function: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
	\item{lastTopVariable}{
		The maximum number of variables to be tested
	}
	\item{timeOutcome}{
		The name of the column in \code{data} that stores the time to event (needed only for a Cox proportional hazards regression model fitting)
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
	\item{z.NeRI}{
		A vector in which each element represents the \emph{z}-score of the NeRI, associated to the \code{testType}, for each feature found in the final model
	}
	\item{loops}{
		The number of loops it took for the model to stabilize
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{updateModel}}}
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
	
	# Rank the variables:
	# - Analyzing the raw data
	# - Using a Cox proportional hazards fitting
	# - According to the NeRI
	rankedDataCancer <- univariateRankVariables(variableList = cancerVarNames,
	                                            formula = "Surv(pgtime, pgstat) ~ 1",
	                                            Outcome = "pgstat",
	                                            data = dataCancer,
	                                            categorizationType = "Raw",
	                                            type = "COX",
	                                            rankingTest = "NeRI",
	                                            description = "Description")
	# Get a Cox proportional hazards model using:
	# - 10 bootstrap loops
	# - The ranked variables
	# - The Wilcoxon rank-sum test as the feature inclusion criterion
	cancerModel <- NeRIBasedFRESA.Model(loops = 10,
	                                    Outcome = "pgstat",
	                                    variableList = rankedDataCancer,
	                                    data = dataCancer,
	                                    type = "COX",
	                                    testType= "Wilcox",
	                                    timeOutcome = "pgtime")
	# Update the model, adding first order interactions
	
	uCancerModel <- updateNeRIModel(Outcome = "pgstat",
	        VarFrequencyTable = cancerModel$ranked.var,
	        variableList = cancerVarNames,
	        data = dataCancer,
	        type = "COX",
	        testType = "Wilcox",
	        timeOutcome = "pgtime",
	        interaction = 2)
	# Shut down the graphics device driver
	dev.off()}
}

\keyword{Model_Generation}