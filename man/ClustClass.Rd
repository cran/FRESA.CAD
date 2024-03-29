\name{ClustClass}
\alias{ClustClass}
\title{Hybrid Hierarchical Modeling}
\description{
	This function returns the outcome associated features and the supervised-classifier present at each one of the unsupervised data clusters
}
\usage{
	ClustClass(formula = formula,
	            data=NULL,
	            filtermethod=univariate_KS,
	            clustermethod=GMVECluster,
	            classmethod=LASSO_1SE,
	            filtermethod.control=list(pvalue=0.1,limit=21),
	            clustermethod.control= list(p.threshold = 0.95, 
	                                   p.samplingthreshold = 0.5),
	            classmethod.control=list(family = "binomial"),
	            pca=TRUE,
	            normalize=TRUE
	            )

}
\arguments{
	\item{formula}{
		An object of class \code{formula} with the formula to be fitted
	}
	\item{data}{
		A data frame where all variables are stored in different columns
	}
	\item{filtermethod}{
		The function name that will return the relevant features
	}
	\item{clustermethod}{
		The function name that will cluster the data points
	}
	\item{classmethod}{
		The function name of the binary classification method
	}
	\item{filtermethod.control}{
		A list with the parameters to be passed to the filter function
	}
	\item{clustermethod.control}{
		A list with the parameters to be passed to the clustering function
	}
	\item{classmethod.control}{
		A list with the parameters to be passed to the classification function
	}
	\item{pca}{
		if TRUE it will compute the PCA transform
	}
	\item{normalize}{
		if pca=TRUE and normalize=TRUE it will normalize all the data.
	}
}
\details{
This function will first call the filter function that should return the relevant a named vector with the p-value of the features associated with the outcome. 
Then it will call user-supplied clustering algorithm that must return a relevant data partition based on the discovered features.
The returned object of the clustering function must contain a $classification object indicates the class label of each data point. 
Finally, the function will call the classification function on each cluster returned by the clustering function.
}
\value{
	\item{features}{
		The named vector of FDR adjusted p-values returned by the filtering function.
	}
	\item{cluster}{
		The clustering function output
	}
	\item{models}{
		The list of classification objects per data cluster
	}
}

\examples{
	\dontrun{	
      library(mlbench) # Location of the Sonar data set
	  library(mclust) # The cluster library
      data(Sonar)
      Sonar$Class <- 1*(Sonar$Class == "M")
	  #Train hierachical classifier
      mc <- ClustClass(Class~.,Sonar,clustermethod=Mclust,clustermethod.control=list(G = 1:4))
	  #report the classification
      pb <- predict(mc,Sonar)
      print(table(1*(pb>0.0),Sonar$Class))
	}
}
\author{Jose G. Tamez-Pena}

\keyword{Model_Generation}