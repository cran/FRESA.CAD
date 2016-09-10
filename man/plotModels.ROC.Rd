\name{plotModels.ROC}
\alias{plotModels.ROC}
\title{Plot test ROC curves of each cross-validation model}
\description{
	This function plots test ROC curves of each model found in the cross validation process.
	It will also aggregate the models into a single prediction performance, plotting the resulting ROC curve (models coherence).
	Furthermore, it will plot the mean sensitivity for a given set of specificities.
}
\usage{
   plotModels.ROC(modelPredictions,
    number.of.models=0,
    specificities=c(0.975,0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),
    theCVfolds=1,
    predictor="Prediction",
	cex=1.0,
    ...)
}
\arguments{
	\item{modelPredictions}{
		A data frame returned by the \code{crossValidationFeatureSelection_Bin} function, either the \code{Models.testPrediction}, the \code{FullBSWiMS.testPrediction},\cr the \code{Models.CVtestPredictions}, the \code{TestRetrained.blindPredictions},\cr the \code{KNN.testPrediction}, or the \code{LASSO.testPredictions} value
	}
	\item{number.of.models}{
		The maximum number of models to plot
	}
	\item{specificities}{
		Vector containing the specificities at which the ROC sensitivities will be calculated
	}
	\item{theCVfolds}{
		The number of folds performed in a Cross-validation experiment
	}
	\item{predictor}{
		The name of the column to be ploted
	}
	\item{cex}{
		Controling the font size of the text inside the plots
	}
	\item{...}{
		Additional parameters for the \code{roc} function (\code{pROC} package)
	}
}
\value{
	\item{ROC.AUCs}{
		A vector with the AUC of each ROC
	}
	\item{mean.sensitivities}{
		A vector with the mean sensitivity at the specificities given by \code{specificities}
	}
	\item{model.sensitivities}{
		A matrix where each row represents the sensitivitiy at the specificities given by \code{specificities} for a different ROC
	}
	\item{specificities}{
		The specificities used to calculate the sensitivities
	}
	\item{senAUC}{
		The AUC of the ROC curve that resulted from using \code{mean.sensitivities}
	}
	\item{predictionTable}{
		The confusion matrix between the outcome and the ensembled prediction
	}
	\item{ensemblePrediction}{
		The ensembled (median prediction) of the repeated predictions
	}
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
	# Rank the variables:
	# - Analyzing the raw data
	# - According to the zIDI
	rankedDataCancer <- univariateRankVariables(variableList = cancerVarNames,
	                                            formula = "Surv(pgtime, pgstat) ~ 1",
	                                            Outcome = "pgstat",
	                                            data = dataCancer,
	                                            categorizationType = "Raw",
	                                            type = "COX",
	                                            rankingTest = "zIDI",
	                                            description = "Description")
	# Get a Cox proportional hazards model using:
	# - The top 7 ranked variables
	# - 10 bootstrap loops in the feature selection procedure
	# - The zIDI as the feature inclusion criterion
	# - 5 bootstrap loops in the backward elimination procedure
	# - A 5-fold cross-validation in the feature selection, 
	#            update, and backward elimination procedures
	# - A 10-fold cross-validation in the model validation procedure
	# - First order interactions in the update procedure
	cancerModel <- crossValidationFeatureSelection_Bin(size = 7,
	                                               loops = 10,
	                                               Outcome = "pgstat",
	                                               timeOutcome = "pgtime",
	                                               variableList = rankedDataCancer,
	                                               data = dataCancer,
	                                               type = "COX",
	                                               selectionType = "zIDI",
	                                               elimination.bootstrap.steps = 5,
	                                               trainRepetition = 5,
	                                               CVfolds = 10,
	                                               interaction = c(1,2))
	# Plot the results of the blind test set predictions made at each 
	# fold of the cross-validation
	cancerModelPlot <- plotModels.ROC(cancerModel$Models.testPrediction)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
