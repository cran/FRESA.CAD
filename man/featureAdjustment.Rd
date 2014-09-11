\name{featureAdjustment}
\alias{featureAdjustment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjust each listed variable to the provided set of covariates.
%%  ~~function to do ... ~~
}
\description{
At each data strata, the adjustment of the variables is done by fitting the candidate variable to the provided linear model on a control population.
If the control residual-variance is reduced significantly, then the residual of the predicted observation to the real observation is used in the resulting data frame.
If the control-residuals are not reduced significantly, then the control-mean is subtracted for all observations.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
featureAdjustment(varList, baseModel, strata = "NA", dataframe, 
referenceframe, type = c("LM", "GLS"), p.value = 0.05, correlationGroup = "ID")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varList}{
  The list of variables to be adjusted
%%     ~~Describe \code{varList} here~~
}
  \item{baseModel}{
  The string that describes the model to be tested at each variable
%%     ~~Describe \code{baseModel} here~~
}
  \item{strata}{
  The column name of the variable with a numeric code of the strata (i.e. 0=Male, 1=Female)
%%     ~~Describe \code{strata} here~~
}
  \item{dataframe}{
  The data frame to be adjusted

}
  \item{referenceframe}{
  The control/reference frame used to fit the model
%%     ~~Describe \code{referenceframe} here~~
}
  \item{type}{
  The type of fitting: LM Linear fitting, or GLS generalized least squares.

}
  \item{p.value}{
  The used p.value to decide if the model will be used to reduce variability
%%     ~~Describe \code{p.value} here~~
}
  \item{correlationGroup}{
  For GLS fitting, the correlationGroup indicates the name of the variable to be used to group the data. 
%%     ~~Describe \code{correlationGroup} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns an adjusted data frame, where each input observation at each strata has been adjusted based on control-fitted model. 

}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Jose G. Tamez-Pena

}
\note{
The function prints the residuals and the F statistic for all candidate variables.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

%%\examples{
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Data_Conditioning}
