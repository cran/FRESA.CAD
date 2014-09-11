\name{backVarNeRIElimination}
\alias{backVarNeRIElimination}
\title{
Residual based backwards elimination 
}
\description{
Remove linear model terms that do not improve the model residuals .
}
\usage{
backVarNeRIElimination(object,
	pvalue=0.05,
	Outcome="Class",
	data,
	startOffset=0, 
	type = c("LOGIT", "LM","COX"),
	testType=c("Binomial","Wilcox","tStudent","Ftest"),
	setIntersect=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  the linear model
}
  \item{pvalue}{
  the maximum term NeRI p.value allowed in the model
}
  \item{Outcome}{
  the Outcome to be predicted by the linear model
}
  \item{data}{
  The data frame that contains the variables
}
  \item{startOffset}{
  Only terms whose model position is larger than the startOffset are candidates to be removed 
}
  \item{type}{
  Fit type: LOGIT: Logistic, LM: Lineal, COX: survival
}
  \item{testType}{
  The test to be performed by the improved residuals function
}
  \item{setIntersect}{
  The model intersect. Set to zero to force a zero intersect.
}


}
\details{
A stage-wise backwards elimination procedure is used to remove model terms that do not improve the residual classification.
Given a model, each term of the model is evaluated according to the change in their residuals. 
For each term, the model residual are compared between the full and reduced model.
The term with the lowest residual improvement is studied at each cycle. If the improvement is not statistically significant 
("Binomial","Wilcox","tStudent" or "Ftest"), then the model term is removed. The procedure is repeated until no terms can be removed.
}
\value{
  \item{back.model }{The final fitted model}
  \item{loops }{Number of loops required}
  \item{reclas.info }{The reclassification data for the final model }
  \item{back.formula }{The character vector of the final model}
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
