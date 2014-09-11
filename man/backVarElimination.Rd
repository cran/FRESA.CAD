\name{backVarElimination}
\alias{backVarElimination}
\title{
IDI-based backwards elimination 
}
\description{
Remove model terms that do not improve the IDI/NRI significantly.
}
\usage{
backVarElimination(object,
	pvalue=0.05,
	Outcome="Class",
	data,
	startOffset=0, 
	type = c("LOGIT", "LM","COX"),
	selectionType=c("zIDI","zNRI"))
}
\arguments{
  \item{object}{ 
  A logistic or COX model
}
  \item{pvalue}{
  the maximum IDI/NRI p.value allowed in the models.  
}
  \item{Outcome}{
  The binary outcome variable
  
}
  \item{data}{
  The data frame that contains all the variables
}
  \item{startOffset}{
  Only terms whose model position is larger than the startOffset are candidates to be removed
  
}
  \item{type}{
  Fit type: LOGIT: Logistic, LM: Lineal, COX: survival
  
}
  \item{selectionType}{
  Selection criteria for feature inclusion: zIDI or zIRI
  
}


}

\details{
For each model term (x_i) the IDI/NRI is computed given the full model. The term with the smallest drop in improvement is selected (largest p-value).
If p(x_i|model)>p.value, then remove the variable. The procedure is repeated until no variable can be removed.
}
\value{
  \item{back.model}{The final model}
  \item{loops}{The number of loops required to remove all non-significant terms}
  \item{reclas.info}{Reports the NRI/IDI statistics of all model terms}
  \item{back.formula}{The character string of the final formula}
  \item{lastRemoved}{Which term was the last to be removed. -1 if non of the terms were removed}
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
\keyword{Model_Generation}
