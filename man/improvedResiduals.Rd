\name{improvedResiduals}
\alias{improvedResiduals}
\title{Estimate the significance of the reduction of predicted residuals}
\description{
	This function will test the hypothesis that, given a set of two residuals (new vs. old), the new ones are better than the old ones as measured with non-parametric tests.
	Four \emph{p}-values are provided: one for the binomial sign test, one for the paired Wilcoxon rank-sum test, one for the paired \emph{t}-test, and one for the \code{F}-test.
	The proportion of subjects that improved their residuals, the proportion that worsen their residuals, and the net residual improvement (NeRI) will be returned.
}
\usage{
	improvedResiduals(oldResiduals,
	                  newResiduals,
	                  testType = c("Binomial", "Wilcox", "tStudent", "Ftest"))
}
\arguments{
	\item{oldResiduals}{
		A vector with the residuals of the original model
	}
	\item{newResiduals}{
		A vector with the residuals of the new model
	}
	\item{testType}{
		Type of non-parametric test to be evaluated: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
}
\details{
	This function will test the hypothesis that the new residuals are "better" than the old residuals.
	To test this hypothesis, four types of tests are performed:
	\enumerate{
		\item The paired \emph{t}-test, which compares the absolute value of the residuals
		\item The paired Wilcoxon rank-sum test, which compares the absolute value of residuals
		\item The binomial sign test, which evaluates whether the number of subjects with improved residuals is greater than the number of subjects with worsened residuals
		\item The \emph{F}-test, which is the standard test for evaluating whether the residual variance is "better" in the new residuals.
	}
	The proportions of subjects that improved and worsen their residuals are returned, and so is the NeRI.
}
\value{
	\item{p1}{
		Proportion of subjects that improved their residuals to the total number of subjects
	}
	\item{p2}{
		Proportion of subjects that worsen their residuals to the total number of subjects
	}
	\item{NeRI}{
		The net residual improvement (\code{p1}-\code{p2})
	}
	\item{p.value}{
		The one tail \emph{p}-value of the test specified in \emph{testType}
	}
	\item{BinP.value}{
		The \emph{p}-value associated with a significant improvement in residuals
	}
	\item{WilcoxP.value}{
		The single sided \emph{p}-value of the Wilcoxon rank-sum test comparing the absolute values of the new and old residuals
	}
	\item{tP.value}{
		The single sided \emph{p}-value of the paired t-test comparing the absolute values of the new and old residuals
	}
	\item{FP.value}{
		The single sided \emph{p}-value of the F-test comparing the residual variances of the new and old residuals
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
	# Get a Cox proportional hazards model using:
	# - 10 bootstrap loops
	# - All variables except for age
	# - The Wilcoxon rank-sum test as the feature inclusion criterion
	cancerModel <- NeRIBasedFRESA.Model(loops = 10,
	                                    Outcome = "pgstat",
	                                    variableList = cancerVarNames[-1,],
	                                    data = dataCancer,
	                                    type = "COX",
	                                    testType= "Wilcox",
	                                    timeOutcome = "pgtime")
	# Add age to the formula of the obtained model
	frm <- format(cancerModel$formula)
	frm[length(frm)] <- paste(frm[length(frm)], "+ age")
	# Fit the new formula to the same data
	cancerModelAge <- modelFitting(formula(frm), dataCancer, "COX")
	# Get the residuals of the original model
	cancerModelRes <- residualForNeRIs(object = cancerModel$final.model,
	                                   testData = dataCancer,
	                                   Outcome = "pgstat")
	# Get the residuals of the model with the added term
	cancerModelAgeRes <- residualForNeRIs(object = cancerModelAge,
	                                      testData = dataCancer,
	                                      Outcome = "pgstat")
	# Estimate the significance of the NeRI when adding age to the model
	NeRI <- improvedResiduals(oldResiduals = cancerModelRes,
	                          newResiduals = cancerModelAgeRes,
	                          testType = "Wilcox")
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Inspection}
