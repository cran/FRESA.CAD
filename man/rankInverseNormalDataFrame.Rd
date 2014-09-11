\name{rankInverseNormalDataFrame}
\alias{rankInverseNormalDataFrame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs a Z transformation of the data using the Rank Inverse Normal Transform
%%  ~~function to do ... ~~
}
\description{
This function takes a data frame and a reference control population to return a z-transformed data set conditioned to the reference population.
Each sample data for each feature column in the data frame is conditionally z-transformed using a inverse normal rank transform.
The rank of the sample in the reference frame is used to get the corresponding z value using the inverse normal transformation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rankInverseNormalDataFrame(varList, dataframe, referenceframe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varList}{
  A data frame that list the variables to be z-transformed
%%     ~~Describe \code{varList} here~~
}
  \item{dataframe}{
  The data frame to be transformed. Where each row has a subject sample and the columns are features.

}
  \item{referenceframe}{
  The reference/control set used to compute the rank of each variable
%%     ~~Describe \code{referenceframe} here~~
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{
Returns a new data frame were each observation has been conditionally z-transformed given the control data
}

%%\references{
%% ~put references to the literature/web site here ~
%%}

\author{
Jose G. Tamez-Pena

}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
## Use the Prostate cancer data
require(rpart)
data(stagec)
dataCancer <- stagec[complete.cases(stagec),]
## List the variables to rank inverse
varlist <- data.frame(c("g2","grade"),c("cells in G2 phase","tumor grade"))
## Set the group of no progression
noProgress <- subset(dataCancer,pgstat==0)
## Rank the variables
CancerZTransform <- rankInverseNormalDataFrame(varlist,dataCancer,noProgress)
## show the heatmap
x <- heatMaps (varlist,outcome='pgstat',dataframe=CancerZTransform,outcomeGain=4)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Data_Conditioning}
