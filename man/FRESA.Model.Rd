\name{FRESA.Model}
\alias{FRESA.Model}
\title{Automated model selection}
\description{
	This function uses a wrapper procedure to select the best features of a non-penalized linear model that best predict the outcome, given the formula of an initial model template (linear, logistic, or Cox proportional hazards), an optimization procedure, and a data frame.
	A filter scheme may be enabled to reduce the search space of the wrapper procedure.
	The false selection rate may be empirically controlled by enabling bootstrapping, and model shrinkage can be evaluated by cross-validation. 
}
\usage{
	FRESA.Model(formula,
	            data,
	            OptType = c("Binary", "Residual"),
	            pvalue = 0.05,
	            filter.p.value = 0.10,
	            loops = 1,
	            maxTrainModelSize = 10,
	            loop.threshold = 20,
	            elimination.bootstrap.steps = 100,
	            bootstrap.steps = 100,
	            interaction = c(1,1),
	            print = TRUE,
	            plots = TRUE,
	            CVfolds = 10,
	            repeats = 1,
	            nk = 0,
	            categorizationType = c("Raw",
	                                   "Categorical",
	                                   "ZCategorical",
	                                   "RawZCategorical",
	                                   "RawTail",
	                                   "RawZTail",
	                                   "Tail"),
	            cateGroups = c(0.1, 0.9),
	            raw.dataFrame = NULL,
	            var.description = NULL,
	            testType = c("zIDI",
	                         "zNRI",
	                         "Binomial",
	                         "Wilcox",
	                         "tStudent",
	                         "Ftest",
	                         "Both"),
	                         zbaggRemoveOutliers=4.0)

}
\arguments{
	\item{formula}{
		An object of class \code{formula} with the formula to be fitted
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{OptType}{
		Optimization type: Based on the integrated discrimination improvement (Binary) index for binary classification ("Binary"), or based on the net residual improvement (NeRI) index for linear regression ("Residual") 
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the \code{testType}, allowed for a term in the model (it will control the false selection rate)
	}
	\item{filter.p.value}{
		The maximum \emph{p}-value, for a variable to be included to the feature selection procedure
	}
	\item{loops}{
		The number of bootstrap loops for the forward selection procedure
	}
	\item{maxTrainModelSize}{
		Maximum number of terms that can be included in the model
	}
	\item{loop.threshold}{
		After \code{loop.threshold} cycles, only variables that have already been selected in previous cycles will be candidates to be selected in posterior cycles
	}
	\item{elimination.bootstrap.steps}{
		The number of bootstrap loops for the backwards elimination procedure
	}
	\item{bootstrap.steps}{
		The number of bootstrap loops for the bootstrap validation procedure
	}
	\item{interaction}{
		A vector of size two. The terms are used by the search and update procedures, respectively. Set to either 1 for first order models, or to 2 for second order models
	}
	\item{print}{
		Logical. If \code{TRUE}, information will be displayed
	}
	\item{plots}{
		Logical. If \code{TRUE}, plots are displayed
	}
	\item{CVfolds}{
		The number of folds for the final cross-validation
	}
	\item{repeats}{
		The number of times that the cross-validation procedure will be repeated
	}
	\item{nk}{
		The number of neighbors used to generate a \emph{k}-nearest neighbors (KNN) classification. If zero, \emph{k} is set to the square root of the number of cases. If less than zero, it will not perform the KNN classification
	}
	\item{categorizationType}{
		How variables will be analyzed: As given in \code{data} ("Raw"); broken into the \emph{p}-value categories given by \code{cateGroups} ("Categorical"); broken into the \emph{p}-value categories given by \code{cateGroups}, and weighted by the \emph{z}-score ("ZCategorical"); broken into the \emph{p}-value categories given by \code{cateGroups}, weighted by the \emph{z}-score, plus the raw values ("RawZCategorical"); raw values, plus the tails ("RawTail"); or raw values, wighted by the \emph{z}-score, plus the tails ("RawZTail")
	}
	\item{cateGroups}{
		A vector of percentiles to be used for the categorization procedure
	}
	\item{raw.dataFrame}{
		A data frame similar to \code{data}, but with unajusted data, used to get the means and variances of the unadjusted data
	}
	\item{var.description}{
		A vector of the same length as the number of columns of \emph{data}, containing a description of the variables
	}
	\item{testType}{
		For an Binary-based optimization, the type of index to be evaluated by the \code{improveProb} function (\code{Hmisc} package): \emph{z}-value of Binary or of NRI. For a NeRI-based optimization, the type of non-parametric test to be evaluated by the \code{improvedResiduals} function: Binomial test ("Binomial"), Wilcoxon rank-sum test ("Wilcox"), Student's \emph{t}-test ("tStudent"), or \emph{F}-test ("Ftest")
	}
	\item{zbaggRemoveOutliers}{
		For linear regresion, zbaggRemoveOutliers is used to set the z-treshold to be used in the outlier detection. 
	}
	
}
\details{
This is the main function of FRESA.CAD given an outcome formula, and a data.frame this function will do an univariate analysis of the data (\code{univariateRankVariables}),
then it will select the top ranked variables; after that it will select the model that best describes the outcome. At output it will return the bootstrapped performance of the model 
(\code{bootstrapValidation_Bin} or \code{bootstrapValidation_Res}). It can be set to report the cross-validation performance of the selection process which will return either 
a  \code{crossValidationFeatureSelection_Bin} or a \code{crossValidationFeatureSelection_Res} object.
}
\value{
	\item{BSWiMS.model}{
		An object of class \code{lm}, \code{glm}, or \code{coxph} containing the final model
	}
	\item{reducedModel}{
		The resulting object of the backward elimination procedure
	}
	\item{univariateAnalysis}{
		A data frame with the results from the univariate analysis
	}
	\item{forwardModel}{
		The resulting object of the feature selection function.
	}
	\item{updatedforwardModel}{
		The resulting object of the the update procedure
	}
	\item{bootstrappedModel}{
		The resulting object of the bootstrap procedure on \code{final.model}
	}
	\item{cvObject}{
		The resulting object of the cross-validation procedure
	}
	\item{used.variables}{
		The number of terms that passed the filter procedure
	}
	\item{call}{
		the function call
	}
	
}
\references{Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: from area under the ROC curve to reclassification and beyond. \emph{Statistics in medicine} \bold{27}(2), 157-172.}
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
	# Get a Cox proportional hazards model using:
	# - 250 bootstrap loops
	md <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
	                  data = dataCancer,
	                  loops = 250,
					  var.description = cancerVarNames[,2])
	# Get a Cox proportional hazards model using:
	# - 250 bootstrap loops
	# - First order interactions in the update procedure
	md <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
	                  data = dataCancer,
	                  loops = 250,
	                  interaction = c(1,2),
					  var.description = cancerVarNames[,2])
	# Get a Cox proportional hazards model using:
	# - No bootstrapping
	# - No cross-validation
	md <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
	                  data = dataCancer,
	                  CVfolds = 0,
	                  elimination.bootstrap.steps = 1,
					  var.description = cancerVarNames[,2])
	# Get a Cox proportional hazards model using:
	# - NeRI-based optimization
	# - 250 bootstrap loops
	# - First order interactions in the update procedure
	md <- FRESA.Model(formula = Surv(pgtime, pgstat) ~ 1,
	                  data = dataCancer,
	                  OptType = "Residual",
	                  loops = 250,
	                  interaction = c(1,2),
					  var.description = cancerVarNames[,2])
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Model_Generation}