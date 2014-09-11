\name{FRESA.Model}
\alias{FRESA.Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Automated model selection
%%  ~~function to do ... ~~
}
\description{
Given an initial model template (Survival or Linear) and optimization procedure (Residual or Classification) and a data frame, this 
function will use a wrapper procedure to select the best features of a non-penalized linear model that best predict the outcome.
A filter scheme may be enabled to reduced the search space of the wrapper procedure. The false selection rate may be empirically controlled 
by enabling bootstrapping, and model shrinkage can be evaluated by cross-validation. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FRESA.Model(formula,
data,
OptType=c("Binary","Residual"),
pvalue=0.05,
filter.p.value=0.10,
loops=1,
maxTrainModelSize=10,
loop.threshold=20,
elimination.bootstrap.steps=100,
bootstrap.steps=100,
interaction=c(1,1),
print=TRUE,
plot=TRUE,
folds=10,
repeats=1,
k=0,
categorizationType=c("Raw","Categorical","ZCategorical",
"RawZCategorical","RawZTail","RawTail"),
cateGroups=c(0.1,0.9),
raw.dataFrame=NULL,
var.description=NULL,
testType=c("zIDI","zNRI","Binomial","Wilcox","tStudent","Ftest")
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
   The Base Formula. It Defines the fitting function
}
  \item{data}{
  The data frame
}
  \item{OptType}{
  The optimization type to be used. \cr
  Binary: Binary classification based on IDI\cr
  Residuals: Residual optimization\cr
}
  \item{pvalue}{
   The maximum p value for FSR-controlled warping procedure 
}
  \item{filter.p.value}{
   The maximum p value for the filter selection procedure.
   }
  \item{loops}{
   The number of bootstrap loops for the forward selection procedure
}
  \item{maxTrainModelSize}{
   The maximum size of the model that can be created
}
  \item{loop.threshold}{
   The loop.threshold is used by the bootstrap procedure to filter-out the non-frequent terms.
}
  \item{elimination.bootstrap.steps}{
   The number of bootstrap steps to be used by the backwards elimination procedure
}
  \item{bootstrap.steps}{
   The number of steps to be performed by the model bootstrapping-validation function.
}
  \item{interaction}{
   The vector of two elements with interactions by the search and update procedures (1: first order mode, 2: Second order models)
}
  \item{print}{
   if TRUE it will print model information
}
  \item{plot}{
	if TRUE it will create the plots
}
  \item{folds}{
	The number of cross-validation folds to be used
}
  \item{repeats}{
  The number of times that the cross-validation procedure will be repeated
}
  \item{k}{
   the number of neighbours for the KNN. If zero k is set to the square root of cases. If KNN is less than 0 it will not perform the KNN
}
  \item{categorizationType}{
  It can be used to create categorical terms on the models
}
  \item{cateGroups}{
  the p-values to be used in for the categorization procedure

}
  \item{raw.dataFrame}{
  The original data frame. It is used to get the mean and variances of the original data
}
  \item{var.description}{
  A vector of the same size as the columns of the data frame that contains a description of the variables.
}
  \item{testType}{
  The test to be used to select a model term 
}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{

%%  If it is a LIST, use
  \item{model}{The final model}
  \item{reducedModel}{The object of the backwards elimination procedure}
  \item{univariateAnalysis}{A Data frame with the univariate analysis}
  \item{firstModel}{The ReclassificationFRESA.Model or NeRIBasedFRESA.Model object }
  \item{updatedModel}{The object after the update procedure}
  \item{bootstrappedModel}{The bootstrapped object of the final model}
  \item{cvObject}{The cross-validation object}
  \item{used.variables}{The number of terms that pass the filter procedure}
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
## Get the stage C prostate cancer data from rpart
library(rpart)
data(stagec)
## Split the stages into several columns
dataCancer <- cbind(stagec[,1:3],stagec[,5:6],1*(stagec[,7] == 4),
  1*(stagec[,7] == 5),1*(stagec[,7] == 6),1*(stagec[,7] == 7),
  1*(stagec[,7] == 8),1*(stagec[,7] >= 9),1*(stagec[,4] == 2),
  1*(stagec[,8]=='diploid'),1*(stagec[,8]=='tetraploid'),1-1*(stagec[,8]=='aneuploid'))
## Assign the column names
colnames(dataCancer) <- c("pgtime","pgstat","age","g2","Grade",
                          "gleason4","gleason5","gleason6","gleason7","gleason8",
                          "gleason910","eet","diploid","tetraploid","Notaneuploid")
## Lets remove the incomplete cases
dataCancer <- as.data.frame(dataCancer[complete.cases(dataCancer),])

## A simple Call for COX modeling
\dontrun{
md <- FRESA.Model(Surv(pgtime,pgstat) ~ 1,dataCancer)
plotModels(md$
}

\dontrun{
## A simple Call for Logistic modeling
md <- FRESA.Model(pgstat ~ 1 ,dataCancer)
}

## A simple Call for LM modeling with residual optimization

\dontrun{
md <- FRESA.Model(pgstat ~ 1 ,dataCancer,OptType="Residual")
}

\dontrun{
## A Call for COX modeling with Bagging
md <- FRESA.Model(Surv(pgtime,pgstat) ~ 1,dataCancer,loops=250)
}

\dontrun{
## A Call for COX modeling with Bagging and second order update interaction
md <- FRESA.Model(Surv(pgtime,pgstat) ~ 1,dataCancer,loops=250,interaction=c(1,2))
}

\dontrun{
## A Call for one COX modeling without Bagging nor Cross validation
md <- FRESA.Model(Surv(pgtime,pgstat) ~ 1,dataCancer,folds=0, 
elimination.bootstrap.steps=1,filter.p.value=0.1)
}

\dontrun{
## A Call for one COX modeling optimizing Matingale residuals with  
## second order interactions in the update procedure
mdr <- FRESA.Model(Surv(pgtime,pgstat) ~ 1,dataCancer,
OptType="Residual",loops=250,interaction=c(1,2))
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Model_Inspection}
