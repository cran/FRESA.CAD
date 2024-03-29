\name{featureAdjustment}
\alias{featureAdjustment}
\title{Adjust each listed variable to the provided set of covariates}
\description{
	This function fits the candidate variables to the provided model formula,for each strata, on a control population.
	If the variance of the residual (the fitted observation minus the real observation) is reduced significantly, then, such residual is used in the resulting data frame.
	Otherwise, the control mean is subtracted to the observation.
}
\usage{
	featureAdjustment(variableList,
	                  baseFormula,
	                  strata = NA,
	                  data,
	                  referenceframe,
	                  type = c("LM", "GLS", "RLM","NZLM","SPLINE","MARS","LOESS"),
	                  pvalue = 0.05,
	                  correlationGroup = "ID",
	                  ...
	                  )
}
\arguments{
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{baseFormula}{
		A string of the type "var1 +...+ varn" that defines the model formula to which variables will be fitted
	}
	\item{strata}{
		The name of the column in \code{data} that stores the variable that will be used to stratify the fitting
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{referenceframe}{
		A data frame similar to \code{data}, but with only the control population
	}
	\item{type}{
		Fit type: linear fitting ("LM"), generalized least squares fitting ("GLS") or Robust ("RLM")
	}
	\item{pvalue}{
		The maximum \emph{p}-value, associated to the \emph{F}-test, for the model to be allowed to reduce variability
	}
	\item{correlationGroup}{
		The name of the column in \code{data} that stores the variable to be used to group the data (only needed if \code{type} defined as "GLS")
	}
	\item{...}{
		parameters for smooth.spline,loess or mda::mars)
	}
}
\value{
	A data frame, where each input observation has been adjusted from \code{data} at each \code{strata}
}
\author{Jose G. Tamez-Pena and Antonio Martinez-Torteya}
\note{This function prints the residuals and the \emph{F}-statistic for all candidate variables}
\keyword{Data_Conditioning}