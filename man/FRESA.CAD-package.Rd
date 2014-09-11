\name{FRESA:CAD-package}
\alias{FRESA:CAD-package}
\alias{FRESA:CAD}
\docType{package}
\title{  
FeatuRE Selection Algorithms for Computer Aided Diagnosis (FRESA.CAD)
}
\description{
FRESA:CAD is a set of tools for the selection and cross-validation of computer-aided-diagnosis/prognosis models based on logistic regression models,
COX models or linear models. The package include functions for univariate analysis, data conditioning (Co-variate adjustment, and normalization), model-selection with false-selection-rate (FSR) control.
and cross-validation of the model selection procedure.
}

\details{
\tabular{ll}{
Package: \tab FRESA:CAD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-28\cr
License: \tab LGPL (>= 2)\cr
}
Purpose: The design of diagnostic or prognostic multivariate models via the selection of significantly discriminant features.
The models are selected via the step-wise selection of features that offer a significant improvement in subject classification/error.
The false selection rate (FSR) is empirical controlled via bootstrapped samples. 
Variables that do not improve subject classification/error on the blind test are not included in the models.

The main function is the selection of a FSR controlled diagnostic/prognostic logistic or COX models constructed from a large set of candidate features.
The variable selection may start by conditioning all variables via covariate-adjustment and z-inverse-rank-transform. 
In order to integrate features with partial discriminant power, the package can be used to categorize the continuous variables and rank its
discriminant power. Once ranked, each feature is bootstrapped-tested in a multivariate model and its blind performance is evaluated.
Variables with statistical significant improvement in classification/error are stored and finally inserted into the final model according to their relative
store frequency. 
A cross-validation procedure may be used to diagnose the amount of model shrinkage produced by the selection scheme.
   
}
\author{
Jose Gerardo Tamez-Pena

Maintainer: <jose.tamezpena@itesm.mx>
}

%%\references{
%%~~ Literature or other references for background information ~~
%%}

\keyword{ package }

%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
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
dataCancer <- dataCancer[complete.cases(dataCancer),]
## Load the list of variables to be searched
data(CancerVarNames)

\dontrun{

## Rank the COX model variables according to the IDI
rankedCancer <- univariateRankVariables(CancerVarNames,"Surv(pgtime,pgstat) ~ 1",
"pgstat", dataCancer,"Raw","COX","zIDI",description="Description");
## Get the bootstrapped model that best fit the prognosis
CancerModel <- ReclassificationFRESA.Model(20,1.0,0.05,10,'1+age','pgstat',rankedCancer,
dataCancer,type="COX",timeOutcome="pgtime")
## Update the model that best fit the prognosis
UCancerModel <- updateModel('pgstat','1+age',0.05,CancerModel$ranked.var,
rankedCancer, dataCancer,type="COX",timeOutcome="pgtime")
## Remove not significant variables
redCancerModel <- backVarElimination(UCancerModel$final.model,pvalue=0.05,
Outcome='pgstat', data=dataCancer,startOffset=0,type="COX");
## Bootstrap the parameters of the final model 
CancerValid <- bootstrapValidation(1.00,50,redCancerModel$back.model$formula,
'pgstat',dataCancer,type="COX");
## Print a summary of the bootstrapped model
s <- summary(CancerValid)
## Plot the model
plot(CancerValid)

## Scale the C prostate cancer data for a heatmap
dataCancersc <- as.data.frame(scale(dataCancer))
## Plot the heat map of all variables
y <- heatMaps (rankedCancer,outcome='pgstat',dataframe=dataCancersc) 
## Plot the heat map of top ranked variables
x <- heatMaps (rankedCancer,CancerModel$ranked.var,outcome='pgstat',
dataframe=dataCancersc,outcomeGain=10) 



}

\dontrun{
## Do a 5-fold cross-validation of a non-bootstrapped model \cr 
modelValidation <- crossValidationFeatureSelection(5,1.0,0.05,1,'1+age',
'pgstat',rankedCancer,dataCancer,type="COX",timeOutcome="pgtime",trainFraction=0.80,
trainRepetition=2,elimination.pValue=0.1)

## List the COX models
modelValidation$formula.list
}
}
