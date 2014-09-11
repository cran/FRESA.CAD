\name{timeSerieAnalysis}
\alias{timeSerieAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the listed time series variables to a given gls model
%%  ~~function to do ... ~~
}
\description{
Plots the time evolution and does a time series analysis of time dependent features columns in a data frame. 
Features listed are fitted to the provided time model (mixed effect model) with the gls procedure.
At output returns for each variable the coefficients, the t values and the corresponding the p.values, 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
timeSerieAnalysis(varList,
baseModel,
dataframe,
timevar="time",
contime=".",
Outcome=".",
...,
description=".",
Ptoshow = c(1),
plegend= c("p"),
timesign= "-",
catgo.names = c("Control", "Case")
) 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varList}{
  The list of variables to be time-modelled with a gls fit

}
  \item{baseModel}{
  the model to be fitted to all the variables in the list
}
  \item{dataframe}{
  The data frame that has column features and time dependent row observations grouped by subject
}
  \item{timevar}{
  The variable that identifies the visit ID of the observations
}
  \item{contime}{
  An optional continuous time variable
}
  \item{Outcome}{
  An optional binary outcome that may be used to show the stratified analysis
}
  \item{...}{
  Addition parameters to be passed to the gls function
}
  \item{description}{
  The name of the column that has the variable description
}
  \item{Ptoshow}{
  List the index of the p values to be shown in the plot
}
  \item{plegend}{
  List the legend of the p values
}
  \item{timesign}{
  The direction of the arrow of time
}
  \item{catgo.names}{
  The legends of the categories binary categories
}

}
\details{
This function will plot the evolution of the mean value of the listed variables with its corresponding error bars.
Then it will fit the data to the provided time model with a generalized least square procedure and it will plot the fitted values. 
If a binary variable was provided, the plots will contain the case and control data. At output the function will return the corresponding
model coefficients, their corresponding t-values, the standard errors and its associated p.values.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{coeff}{The matrix of the gls coefficients'}
  \item{std.Errors}{The matrix of the standardized error of each coefficient}
  \item{t.values}{The matrix of the t.values}
  \item{p.values}{The p.values of the coefficients}
  \item{sigmas}{The RMS value of the fitting}
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

%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Data_Inspection}
