\name{reportEquivalentVariables}
\alias{reportEquivalentVariables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Report the set of variables that will perform an equivalent IDI discriminant function
%%  ~~function to do ... ~~
}
\description{
Given a model, the function will report in a matrix of all variables that may be used in the linear model without affecting its classification performance. 
For each variable in the model, this function will loop all candidate variables and report
all of them that will have an equivalent or better zIDI than the original model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reportEquivalentVariables(object,
pvalue=0.05,
dataframe,
searchVarList,
Outcome="Class", 
type = c("LOGIT", "LM","COX"),
eqFrac=0.9,
description=".") 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  a formula for the fitting function
}
  \item{pvalue}{
  the maximum pvalue for variable equivalence
}
  \item{dataframe}{
  the dataframe to explore
}
  \item{searchVarList}{
  the list of candidate variables
}
  \item{Outcome}{
  the binary outcome
}
  \item{type}{
  the fitting type
}
  \item{eqFrac}{
  a factor to relax the current z-value
}
  \item{description}{
  The column with the variables description
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{
Returns a data frame with three columns. The first column is the original variable of the model. 
The second column list all variables that will not statistically affect the original model performance 
The third column list the corresponding zIDI for each equivalent variable.
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
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Inspection}
