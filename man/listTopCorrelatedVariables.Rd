\name{listTopCorrelatedVariables}
\alias{listTopCorrelatedVariables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List the variables that are highly correlated
%%  ~~function to do ... ~~
}
\description{
Computes the correlation (Pearson, Spearman or Kendall) for each listed variable in the data set and returns for each variable a list of the
variables that are correlated. It also provides a short variable list without the higly correlated variables.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
listTopCorrelatedVariables(variablesNames,dataframe,
					pvalue=0.001,
					corthreshold=0.9,
					method = c("pearson", "kendall", "spearman"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variablesNames}{ 
  the data frame with the list of variables to be explored
}
  \item{dataframe}{
  The data frame that contains the variables
}
  \item{pvalue}{
  the maximum p.value to be used to call a pair of variables significantly correlated
}
  \item{corthreshold}{
  The minimum correlation threshold
}  \item{method}{
  The correlation method: Pearson, Spearman or Kendall
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{

	\item{correlated.variables}{ Dataframe with two columns: 1) cor.var.names: The variables that are correlated 2) cor.var.value: The correlation value}
	\item{short.list}{ A vector with the short list of variables that first enter the correlation analysis. Correlated variables that were correlated to the short list are not included.}

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
