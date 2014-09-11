\name{heatMaps}
\alias{heatMaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the heat map of the selected variables
%%  ~~function to do ... ~~
}
\description{
Create a heat map for a data set based on the univariate ranking or frequency ranking
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heatMaps(varList,
varRank=NULL,
outcome,
idnames=NULL,
dataframe,
title="Heat Map",
colCluster="NA",
prediction=NULL,
outcomeGain=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varList}{
   A data frame with a list of variables to be used in the heat map
}
  \item{varRank}{
  a data frame with the rank of the variables
}
  \item{outcome}{
   the column name of the outcome variable to be inserted into the heat map 
}
  \item{idnames}{
   the column with the ids of the measurements
}
  \item{dataframe}{
   the dataframe to be used to extract the heatmap
}
  \item{title}{
   the title of the heat map
}
  \item{colCluster}{
   Set to true to cluster the variables
}
  \item{prediction}{
   A vector with a prediction for each subjects. This will be used to rank the heatmap
}
  \item{outcomeGain}{
   An optional constant that can be used to force the heatmap cluster subjects according to the outcome.
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

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

\value{

%%  If it is a LIST, use
  \item{dataMatrix}{A data matrix of all the listed terms described by varList }
  \item{orderMatrix}{The data matrix ordered according to the outcome}
  \item{heatMap}{The heat map object}
  
 }

%%\examples{
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Data_Inspection}
