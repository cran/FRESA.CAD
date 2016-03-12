\name{ForwardSelection.Model.Res}
\alias{ForwardSelection.Model.Res}
\title{NeRI-based feature selection procedure for linear, logistic, or Cox proportional hazards regression models}
\description{
	This function performs a bootstrap sampling to rank the most frequent variables that statistically aid the models by minimizing the residuals.
	After the frequency rank, the function uses a forward selection procedure to create a final model, whose terms all have a significant contribution to the net residual improvement (NeRI).
}
\usage{
	ForwardSelection.Model.Res(size = 100, 
	                     fraction = 1, 
	                     pvalue = 0.05, 
	                     loops = 100, 
	                     covariates = "1", 
	                     Outcome, 
	                     variableList, 
	                     data, 
	                     maxTrainModelSize = 10, 
	                     type = c("LM", "LOGIT", "COX"), 
	                     testType=c("Binomial", "Wilcox", "tStudent", "Ftest"),
	                     timeOutcome = "Time",
	                     loop.threshold = 20,
	                     interaction = 1,
	                     cores = 4)
}
\arguments{
	\item{size}{
		The number of candidate variables to be tested (the first \code{size} variables from \code{variableList})
	}
	\item{fraction}{
		The fraction of data (sampled with replacement) to be used as train
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the NeRI, allowed for a term in the model (controls the false selection rate)
	}
	\item{loops}{
		The number of bootstrap loops
	}
	\item{covariates}{
		A string of the type "1 + var1 + var2" that defines which variables will always be included in the models (as covariates)
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{maxTrainModelSize}{
		Maximum number of terms that can be included in the model
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{testType}{
		Type of non-parametric test to be evaluated by the \code{improvedResiduals} function: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
	\item{timeOutcome}{
		The name of the column in \code{data} that stores the time to event (needed only for a Cox proportional hazards regression model fitting)
	}
	\item{loop.threshold}{
		After \code{loop.threshold} cycles, only variables that have already been selected in previous cycles will be candidates to be selected in posterior cycles
	}
	\item{interaction}{
		Set to either 1 for first order models, or to 2 for second order models
	}
	\item{cores}{
		Cores to be used for parallel processing
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
	\item{ranked.var}{
		An array with the ranked frequencies of the features
	}
	\item{z.NeRIs}{
		A vector in which each element represents the \emph{z}-score of the NeRI, associated to the \code{testType}, for each feature found in the final model
	}
	\item{formula.list}{
		A list containing objects of class \code{formula} with the formulas used to fit the models found at each cycle
	}
	\item{variableList}{
		A list of variables used in the forward selection
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{ForwardSelection.Model.Bin}}}
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
	cancerModel <- ForwardSelection.Model.Res(loops = 10,
	                                    Outcome = "pgstat",
	                                    variableList = rankedDataCancer,
	                                    data = dataCancer,
	                                    type = "COX",
	                                    testType= "Wilcox",
	                                    timeOutcome = "pgtime")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Generation}
