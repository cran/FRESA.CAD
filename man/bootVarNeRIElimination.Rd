\name{bootVarNeRIElimination}
\alias{bootVarNeRIElimination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove model terms that do not improve blind-test residuals
%%  ~~function to do ... ~~
}
\description{
Remove model terms that do not improve test-prediction residuals based on a bootstrap procedure. 
Test Residuals are evaluated by using the bootstrap validation procedure. 
The non-observed set of the bootstrap sample is used to estimate the test-residuals on the full and
the reduced model. Terms that with the smallest improvement in predicted test residuals may 
be removed if the p-value("Binomial","Wilcox","tStudent" or "Ftest") is greater than the specified pvalue
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootVarNeRIElimination(object,
	pvalue=0.05,
	Outcome="Class",
	data,
	startOffset=0, 
	type = c("LOGIT", "LM","COX"),
	testType=c("Binomial","Wilcox","tStudent","Ftest"),
	bootLoops=250,
	bootFraction=1.0,
	setIntersect=1
	)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  the linear model
}
  \item{pvalue}{
  the critical p.value
  
}
  \item{Outcome}{
  the name of the outcome
  
}
  \item{data}{
  The data frame that contains the variables
}
  \item{startOffset}{
  Start removing from position
  
}
  \item{type}{
  Fit type: LOGIT: Logistic, LM: Lineal, COX: survival
  
}
  \item{testType}{
  Test to be used to compare the residuals of the full and reduced model: "Binomial","Wilcox","tStudent","Ftest"

}
  \item{bootLoops}{
  The number of bootstrap samples

}
  \item{bootFraction}{
  the sampled fraction used for the train-independent random data sampling

}
  \item{setIntersect}{
  The intersect of the model

}

}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{

%%  If it is a LIST, use
  \item{back.model }{The model}
  \item{loops }{The number of loops performed}
  \item{bootCV }{The bootstrap of the formula}
  \item{back.formula }{The final formula}
  \item{lastRemoved }{The last removed item}
%% ...
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
%%
%%
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Generation}
