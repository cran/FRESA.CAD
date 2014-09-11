\name{bootstrapVarElimination}
\alias{bootstrapVarElimination}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Backwards elimination of no-significant model terms
%%  ~~function to do ... ~~
}
\description{
Removes model terms that do improve to the subject classification. 
The improvement in subject classification for each term is estimated by the bootstrapValidation procedure.
The blind set of each bootstrap sample is used to estimate the improvement in subject classification.
The term with the smallest improvement is selected and If that term does not 
significantly improve the bootstrapped subject classification, then the term is removed.   
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrapVarElimination(object,
	pvalue=0.05,
	Outcome="Class",
	data,
	startOffset=0, 
	type = c("LOGIT", "LM","COX"),
	selectionType=c("zIDI","zNRI"),
	bootLoops=250,
	bootFraction=1.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  the linear model
}
  \item{pvalue}{
  the desired maximum p.value
  
}
  \item{Outcome}{
  the desired maximum p.value
  
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
  \item{selectionType}{
  Selection criteria for feature inclusion: zIDI or zIRI

}
  \item{bootLoops}{
  The number of replicas for the bootstrap estimation of parameters

}
  \item{bootFraction}{
  the sampled fraction used for the train-independent random data sampling

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
