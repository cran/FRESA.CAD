\name{bootstrapVarElimination_Bin}
\alias{bootstrapVarElimination_Bin}
\title{IDI/NRI-based backwards variable elimination with bootstrapping}
\description{
	This function removes model terms that do not improve the bootstrapped integrated discrimination improvement (IDI) or net reclassification improvement (NRI) significantly.
}
\usage{
	bootstrapVarElimination_Bin(object,
	                        pvalue = 0.05,
	                        Outcome = "Class",
	                        data,
	                        startOffset = 0, 
	                        type = c("LOGIT", "LM", "COX"),
	                        selectionType = c("zIDI", "zNRI"),
	                        loops = 250,
	                        fraction = 1.0,
	                        print=TRUE,
	                        plots=TRUE,
							adjsize=1)
}
\arguments{
	\item{object}{ 
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the model to be analyzed
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to either IDI or NRI, allowed for a term in the model
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
	\item{selectionType}{
		The type of index to be evaluated by the \code{improveProb} function (\code{Hmisc} package): \emph{z}-score of IDI or of NRI
	}
	\item{loops}{
		The number of bootstrap loops
	}
	\item{fraction}{
		The fraction of data (sampled with replacement) to be used as train
	}
	\item{print}{
		Logical. If \code{TRUE}, information will be displayed
	}
	\item{plots}{
		Logical. If \code{TRUE},  plots are displayed
	}
	\item{adjsize}{
		the number of features to be used in the BH FDR correction
	}
}
\details{
	For each model term \eqn{x_i}, the IDI or NRI is computed for the Full model and the reduced model( where the term \eqn{x_i} removed).
	The term whose removal results in the smallest drop in bootstrapped improvement is selected. The hypothesis: the 
	term adds classification improvement is tested by checking the pvalue of average improvement. If \eqn{p(IDI or NRI)>pvalue}, then the term is removed. 
	In other words, only model terms that significantly aid in subject classification are kept.
	The procedure is repeated until no term fulfils the removal criterion.
}
\value{
	\item{back.model}{
		An object of the same class as \code{object} containing the reduced model
	}
	\item{loops}{
		The number of loops it took for the model to stabilize
	}
	\item{reclas.info}{
		A list with the NRI and IDI statistics of the reduced model, as given by the \code{getVar.Bin} function
	}
	\item{bootCV}{
		An object of class \code{bootstrapValidation_Bin} containing the results of the bootstrap validation in the reduced model
	}
	\item{back.formula}{
		An object of class \code{formula} with the formula used to fit the reduced model
	}
	\item{lastRemoved}{
		The name of the last term that was removed (-1 if all terms were removed)
	}
	\item{beforeFSC.model}{
		the beforeFSC model will have the model with the minimum bootstrap test error
	}
	\item{beforeFSC.formula}{
		the string formula of the model used to find the minimum bootstrap test error
	}
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\seealso{\code{\link{bootstrapVarElimination_Res},
				 \link{backVarElimination_Bin},
				 \link{backVarElimination_Res}}}
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
	# - zIDI as the feature inclusion criterion
	# - First order interactions
	cancerModel <- ForwardSelection.Model.Bin(pvalue = 0.1,
	                                           loops = 10,
	                                           covariates = "1 + age",
	                                           Outcome = "pgstat",
	                                           variableList = cancerVarNames,
	                                           data = dataCancer,
	                                           type = "COX",
	                                           timeOutcome = "pgtime",
	                                           selectionType = "zIDI",
	                                           interaction = 2)
	# Remove not significant variables from the previous model:
	# - Using a strict p-value
	# - Excluding the covariate as a candidate for feature removal
	# - Using zIDI as the feature removal criterion
	# - Using 50 bootstrap loops
	reducedCancerModel <- bootstrapVarElimination_Bin(object = cancerModel$final.model,
	                                              pvalue = 0.005,
	                                              Outcome = "pgstat",
	                                              data = dataCancer,
	                                              startOffset = 1,
	                                              type = "COX",
	                                              selectionType = "zIDI",
	                                              loops = 50)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Generation}