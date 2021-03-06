\name{nearestCentroid}
\alias{nearestCentroid}
\title{Class Label Based on the Minimum Mahalanobis Distance}
\description{
   The function will return the set of labels of a data set
}
\usage{
   nearestCentroid(dataset,
                  clustermean=NULL,
                  clustercov=NULL, 
                  p.threshold=1.0e-6)

}
\arguments{
   \item{dataset}{
      The data set to be labeled
   }
   \item{clustermean}{
      The list of cluster centers.
   }
   \item{clustercov}{
      The list of cluster covariances
   }
   \item{p.threshold}{
      The minimum aceptance p.value
   }
}
\details{ The data set will be labeled based on the nearest cluster label. Points distance with membership probability
lower than the acceptance threshold will have the "0" label.
}
\value{
   \item{ClusterLabels}{
      The labels of each point
   }
}

\author{Jose G. Tamez-Pena}

\keyword{Data Labeling}