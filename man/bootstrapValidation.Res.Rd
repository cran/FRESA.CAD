\name{bootstrapValidation_Res}
\alias{bootstrapValidation_Res}
\title{Bootstrap validation of regression models}
\description{
	This function bootstraps the model \emph{n} times to estimate for each variable the empirical bootstrapped distribution of model coefficients, and net residual improvement (NeRI).
	At each bootstrap the non-observed data is predicted by the trained model, and statistics of the test prediction are stores and reported.
}
\usage{
	bootstrapValidation_Res(fraction = 1,
	                        loops = 200,
	                        model.formula,
	                        Outcome,
	                        data,
	                        type = c("LM", "LOGIT", "COX"),
	                        plots = TRUE)
}
\arguments{
	\item{fraction}{
		The fraction of data (sampled with replacement) to be used as train
	}
	\item{loops}{
		The number of bootstrap loops
	}
	\item{model.formula}{
		An object of class \code{formula} with the formula to be used
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{plots}{
		Logical. If \code{TRUE}, density distribution plots are displayed
	}
}
\details{
	The bootstrap validation will estimate the confidence interval of the model coefficients and the NeRI.
	It will also compute the train and blind test root-mean-square error (RMSE), as well as the distribution of the NeRI \emph{p}-values.
}
\value{
	\item{data}{
		The data frame used to bootstrap and validate the model
	}
	\item{outcome}{
		A vector with the predictions made by the model
	}
	\item{boot.model}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing a model whose coefficients are the median of the coefficients of the bootstrapped models
	}
	\item{NeRIs}{
		A matrix with the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{tStudent.pvalues}{
		A matrix with the \emph{t}-test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap train sets
	}
	\item{wilcox.pvalues}{
		A matrix with the Wilcoxon rank-sum test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap train sets
	}
	\item{bin.pvlaues}{
		A matrix with the binomial test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap train sets
	}
	\item{F.pvlaues}{
		A matrix with the \emph{F}-test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap train sets
	}
	\item{test.tStudent.pvalues}{
		A matrix with the \emph{t}-test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{test.wilcox.pvalues}{
		A matrix with the Wilcoxon rank-sum test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{test.bin.pvlaues}{
		A matrix with the binomial test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{test.F.pvlaues}{
		A matrix with the \emph{F}-test \emph{p}-value of the NeRI for each model term, estimated using the bootstrap test sets
	}
	\item{testPrediction}{
		A vector that contains all the individual predictions used to validate the model in the bootstrap test sets
	}
	\item{testOutcome}{
		A vector that contains all the individual outcomes used to validate the model in the bootstrap test sets
	}
	\item{testResiduals}{
		A vector that contains all the residuals used to validate the model in the bootstrap test sets
	}
	\item{trainPrediction}{
		A vector that contains all the individual predictions used to validate the model in the bootstrap train sets
	}
	\item{trainOutcome}{
		A vector that contains all the individual outcomes used to validate the model in the bootstrap train sets
	}
	\item{trainResiduals}{
		A vector that contains all the residuals used to validate the model in the bootstrap train sets
	}
	\item{testRMSE}{
		The global RMSE, estimated using the bootstrap test sets
	}
	\item{trainRMSE}{
		The global RMSE, estimated using the bootstrap train sets
	}
	\item{trainSampleRMSE}{
		A vector with the RMSEs in the bootstrap train sets
	}
	\item{testSampledRMSE}{
		A vector with the RMSEs in the bootstrap test sets
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{bootstrapValidation_Bin},
				 \link{plot.bootstrapValidation_Res}}}
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
	# - Age as a covariate
	# - The Wilcoxon rank-sum test as the feature inclusion criterion
	cancerModel <- ForwardSelection.Model.Res(loops = 10,
	                                    covariates = "1 + age",
	                                    Outcome = "pgstat",
	                                    variableList = cancerVarNames,
	                                    data = dataCancer,
	                                    type = "COX",
	                                    testType= "Wilcox",
	                                    timeOutcome = "pgtime")
	# Validate the previous model:
	# - Using 50 bootstrap loops
	bootCancerModel <- bootstrapValidation_Res(loops = 50,
	                                           model.formula = cancerModel$formula,
	                                           Outcome = "pgstat",
	                                           data = dataCancer,
	                                           type = "COX")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Diagnosis}