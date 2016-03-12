\name{baggedModel}
\alias{baggedModel}
\title{Get the bagged model from a list of forward models}
\description{
	This function will take the frequency-ranked of variables and the list of models to create a single bagged model
}
\usage{
	baggedModel(modelFormulas,
				data,
				type=c("LM","LOGIT","COX"),
				Outcome=NULL,
				timeOutcome=NULL,
				pvalue=0.05,
				backElimination=FALSE,
				frequencyThreshold=0.05,
				removeOutliers=4.0
				)
}
\arguments{
	\item{modelFormulas}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{data}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{type}{
		Fit type: Logistic ("LOGIT"), linear ("LM"), or Cox proportional hazards ("COX")
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the time to outcome
	}
	\item{timeOutcome}{
		The name of the column in \code{data} that stores the time to event (needed only for a Cox proportional hazards regression model fitting)
	}
	\item{pvalue}{
		The elimination p-value)
	}
	\item{backElimination}{
		set it to TRUE if backelimination will be performed at each formula before bagging the coefficients)
	}
	\item{frequencyThreshold}{
		set the frequency the thresold of the frequence of features to be included in the model)
	}
	\item{removeOutliers}{
		The z value for removing outliers from data set)
	}
}
\value{
  \item{bagged.model}{
		the bagged model
	}
	\item{formula}{
		the formula of the model
	}
	\item{frequencyTable}{
		the table of variables ranked by their model frequency
	}
	\item{faverageSize}{
		the average size of the models
	}
	\item{zvalues}{
		The average z-values of the model coefficients
	}
	\item{reducedDataSet}{
		A data set with the outliers removed
	}
	\item{MAD}{
		The mean absolute difference(MAD) of the residuals
	}
	
}
\author{Jose G. Tamez-Pena}
\seealso{\code{\link{medianPredict}}}

\keyword{Model_Generation}