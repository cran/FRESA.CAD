\name{heatMaps}
\alias{heatMaps}
\title{Plot a heat map of selected variables}
\description{
	This function creates a heat map for a data set based on a univariate or frequency ranking
}
\usage{
	heatMaps(variableList,
	         varRank = NULL,
	         Outcome,
	         idnames = NULL,
	         data,
	         title = "Heat Map",
	         colCluster = FALSE,
	         prediction = NULL,
	         outcomeGain = 0)
}
\arguments{
	\item{variableList}{
		A data frame with two columns. The first one must have the names of the candidate variables and the other one the description of such variables
	}
	\item{varRank}{
		A data frame with the name of the variables in \code{variableList}, ranked according to a certain metric
	}
	\item{Outcome}{
		The name of the column in \code{data} that stores the variable to be predicted by the model
	}
	\item{idnames}{
		The name of the column in \code{data} that stores the IDs of the measurements (optional)
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{title}{
		The title of the plot
	}
	\item{colCluster}{
		Logical. If \code{TRUE}, variables will be clustered
	}
	\item{prediction}{
		A vector with a prediction for each subject, which will be used to rank the heat map
	}
	\item{outcomeGain}{
		An optional value that can be used to force the heatmap to cluster subjects according to \code{outcome}
	}
}
\value{
	\item{dataMatrix}{
		A matrix with all the terms in \code{data} described by \code{variableList}
	}
	\item{orderMatrix}{
		A matrix similar to \code{dataMatrix}, where rows are ordered according to the \code{outcome}
	}
	\item{heatMap}{
		A list with the values returned by the \code{heatmap.2} function (\code{gplots} package)
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
	# Scale the C prostate cancer data for a heatmap
	dataCancerScale <- as.data.frame(scale(dataCancer))
	# Generate a heat map using:
	# - The top ranked variables
	# - The scaled data
	hmTop <- heatMaps(variableList = rankedDataCancer,
	                  varRank = cancerModel$ranked.var,
	                  Outcome = "pgstat",
	                  data = dataCancerScale,
	                  outcomeGain = 10)
	# Shut down the graphics device driver
	dev.off()}
}
\keyword{Data_Inspection}
