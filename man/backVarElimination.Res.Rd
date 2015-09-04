\name{backVarElimination_Res}
\alias{backVarElimination_Res}
\title{NeRI-based backwards variable elimination}
\description{
	This function removes model terms that do not significantly improve the "net residual" (NeRI) 
}
\usage{
	backVarElimination_Res(object,
	                       pvalue = 0.05,
	                       Outcome = "Class",
	                       data,
	                       startOffset = 0, 
	                       type = c("LOGIT", "LM", "COX"),
	                       testType = c("Binomial", "Wilcox", "tStudent", "Ftest"),
	                       setIntersect = 1,
						   adjsize= 1)
}
\arguments{
	\item{object}{ 
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the NeRI, allowed for a term in the model
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{startOffset}{
		Only terms whose position in the model is larger than the \code{startOffset} are candidates to be removed
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{testType}{
		Type of non-parametric test to be evaluated by the \code{improvedResiduals} function: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
	\item{setIntersect}{
		The intersect of the model (To force a zero intersect, set this value to 0)
	}
	\item{adjsize}{
		The number of features to be used in the BH FSR correction
	}
}
\details{
	For each model term \eqn{x_i}, the residuals are computed for the Full model and the reduced model( where the term \eqn{x_i} removed).
	The term whose removal results in the smallest drop in residuals improvement is selected. The hypothesis: the 
	term improves residuals is tested by checking the pvalue of improvement. If \eqn{p(residuals better than reduced residuals)>pvalue}, then the term is removed. 
	In other words, only model terms that significantly aid in improving residuals are kept.
	The procedure is repeated until no term fulfils the removal criterion.
	The p-values of improvement can be computed via a sign-test (Binomial) a paired Wilcoxon test, paired t-test or f-test. The first three tests compare the absolute values of
	the residuals, while the f-test test if the variance of the residuals is improved significantly.
}
\value{
	\item{back.model}{
		An object of the same class as \code{object} containing the reduced model
	}
	\item{loops}{
		The number of loops it took for the model to stabilize
	}
	\item{reclas.info}{
		A list with the NeRI statistics of the reduced model, as given by the \code{getVar.Res} function
	}
	\item{back.formula}{
		An object of class \code{formula} with the formula used to fit the reduced model
	}
	\item{lastRemoved}{
		The name of the last term that was removed (-1 if all terms were removed)
	}
	\item{beforeFSC.model}{
		the model with before the FSR procedure. Coefficients are bagged 
	}
	\item{beforeFSC.formula}{
		the string formula of the the FSR procedure
	}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{backVarElimination_Bin},
				 \link{bootstrapVarElimination_Bin}
				 \link{bootstrapVarElimination_Res}}}
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
	# - A lax p-value
	# - 10 bootstrap loops
	# - Age as a covariate
	# - The Wilcoxon rank-sum test as the feature inclusion criterion
	cancerModel <- ForwardSelection.Model.Res(pvalue = 0.1,
	                                    loops = 10,
	                                    covariates = "1 + age",
	                                    Outcome = "pgstat",
	                                    variableList = cancerVarNames,
	                                    data = dataCancer,
	                                    type = "COX",
	                                    testType= "Wilcox",
	                                    timeOutcome = "pgtime")
	# Remove not significant variables from the previous model:
	# - Using a strict p-value
	# - Excluding the covariate as a candidate for feature removal 
	# - Using the Wilcoxon rank-sum test as the feature removal criterion
	reducedCancerModel <- backVarElimination_Res(object = cancerModel$final.model,
	                                             pvalue = 0.005,
	                                             Outcome = "pgstat",
	                                             data = dataCancer,
	                                             startOffset = 1,
	                                             type = "COX",
	                                             testType = "Wilcox")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Generation}