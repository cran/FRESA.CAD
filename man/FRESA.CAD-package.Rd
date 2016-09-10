\name{FRESA.CAD-package}
\alias{FRESA.CAD-package}
\alias{FRESA.CAD}
\docType{package}
\title{FeatuRE Selection Algorithms for Computer-Aided Diagnosis (FRESA.CAD)}
\description{
	Contains a set of utilities for building and testing formula-based models for Computer Aided Diagnosis/prognosis applications via feature selection.
	Bootstrapped Stage Wise Model Selection (B:SWiMS) controls the false selection (FS) for linear, logistic, or Cox proportional hazards regression models.
	Utilities include functions for: univariate/longitudinal analysis, data conditioning (i.e. covariate adjustment and normalization), model validation and visualization.
}	


\details{
	\tabular{ll}{
		Package: \tab FRESA.CAD\cr
		Type: \tab Package\cr
		Version: \tab 2.2.1\cr
		Date: \tab 2016-4-18\cr
		License: \tab LGPL (>= 2)\cr
	}
	Purpose: The design of diagnostic or prognostic multivariate models via the selection of significantly discriminant features.
	The models are selected via the bootstrapped step-wise selection of model features that offer a significant improvement in subject classification/error.
	The false selection control is achieved by train-test partitions, where train sets are used to select variables and test sets used to evaluate model performance.
	Variables that do not improve subject classification/error on the blind test are not included in the models.

	The main function of this package is the selection and cross-validation of diagnostic/prognostic linear, logistic, or Cox proportional hazards regression model constructed from a large set of candidate features.
	The variable selection may start by conditioning all variables via a covariate-adjustment and a \emph{z}-inverse-rank-transformation. 
	In order to integrate features with partial discriminant power, the package can be used to categorize the continuous variables and rank their discriminant power.
	Once ranked, each feature is bootstrap-tested in a multivariate model, and its blind performance is evaluated.
	Variables with a statistical significant improvement in classification/error are stored and finally inserted into the final model according to their relative store frequency. 
	A cross-validation procedure may be used to diagnose the amount of model shrinkage produced by the selection scheme.
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{
Jose Gerardo Tamez-Pena, Antonio Martinez-Torteya and Israel Alanis\cr
Maintainer: <jose.tamezpena@itesm.mx>
}
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
	# - The default parameters
	md <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
	                  data = dataCancer,
					  var.description = cancerVarNames[,2])
	# Get a logistic regression model using
	# - The default parameters
	md <- FRESA.Model(formula = pgstat ~ 1,
	                  data = dataCancer,
					  var.description = cancerVarNames[,2])
	# Get a logistic regression model using:
	# - redidual-based optimization
	md <- FRESA.Model(formula = pgstat ~ 1,
	                  data = dataCancer,
	                  OptType = "Residual",
					  var.description = cancerVarNames[,2])
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
	# - 10 bootstrap loops
	# - Age as a covariate
	# - zIDI as the feature inclusion criterion
	cancerModel <- ForwardSelection.Model.Bin(loops = 10,
	                                           covariates = "1 + age",
	                                           Outcome = "pgstat",
	                                           variableList = rankedDataCancer,
	                                           data = dataCancer,
	                                           type = "COX",
	                                           timeOutcome = "pgtime",
	                                           selectionType = "zIDI")
	# Update the model
	uCancerModel <- updateModel.Bin(Outcome = "pgstat",
	                            VarFrequencyTable = cancerModel$ranked.var,
	                            variableList = rankedDataCancer,
	                            data = dataCancer,
	                            type = "COX",
	                            timeOutcome = "pgtime")
	# Remove not significant variables from the previous model:
	# - Using zIDI as the feature removal criterion
	reducedCancerModel <- backVarElimination_Bin(object = uCancerModel$final.model,
	                                         Outcome = "pgstat",
	                                         data = dataCancer,
	                                         type = "COX",
	                                         selectionType = "zIDI")
	# Validate the previous model:
	# - Using 50 bootstrap loops
	bootCancerModel <- bootstrapValidation_Bin(loops = 50,
	                                       model.formula = reducedCancerModel$back.formula,
	                                       Outcome = "pgstat",
	                                       data = dataCancer,
	                                       type = "COX")	
	# Get the summary of the bootstrapped model
	sumBootCancerModel <- summary.bootstrapValidation_Bin(object = bootCancerModel)
	# Plot the bootstrap results
	plot(bootCancerModel)
	# Scale the C prostate cancer data
	dataCancerScale <- as.data.frame(scale(dataCancer))
	# Generate a heat map using:
	# - All the variables
	# - The scaled data
	hmAll <- heatMaps(variableList = rankedDataCancer,
	                  Outcome = "pgstat",
	                  data = dataCancerScale,
	                  Scale = 10)
	# Generate a heat map using:
	# - The top ranked variables
	# - The scaled data
	hmTop <- heatMaps(variableList = rankedDataCancer,
	                  varRank = cancerModel$ranked.var,
	                  Outcome = "pgstat",
	                  data = dataCancerScale,
	                  Scale = 10)
	# Get a new Cox proportional hazards model using:
	# - The top 5 ranked variables
	# - No bootstrapping
	# - Age as a covariate
	# - The zIDI as the feature inclusion criterion
	# - A train fraction of 0.8
	# - A 2-fold cross-validation in the feature selection and update procedures
	# - A 10-fold cross-validation in the model validation procedure
	# - An elimination p-value of 0.1
	cancerModelCV <- crossValidationFeatureSelection_Bin(size = 5,
	                                                 loops = 1,
	                                                 covariates = "1 + age",
	                                                 Outcome = "pgstat",
	                                                 timeOutcome = "pgtime",
	                                                 variableList = rankedDataCancer,
	                                                 data = dataCancer,
	                                                 type = "COX",
	                                                 selectionType = "zIDI",
	                                                 trainFraction = 0.8,
	                                                 trainRepetition = 2,
	                                                 CVfolds = 10,
	                                                 elimination.pValue = 0.1)
	# List the COX models
	cancerModelCV$formula.list
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{package}