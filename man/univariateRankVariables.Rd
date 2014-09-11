\name{univariateRankVariables}
\alias{univariateRankVariables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate analysis of data 
%%  ~~function to do ... ~~
}
\description{
The function reports the mean, the standard deviation for each listed variable data.
It does a KS test on the raw and z-standardized data. After that, reports the raw t test, the z-standardized t-test, 
the Willcox test, the NRI, the IDI, NeRI and the AUC. Furthermore, it reports the z-value of the variable significance
on the regression model fitted either by LOGIT, LM or COX procedures. The results of the analyses is sorted according the specified test.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
univariateRankVariables(varList, baseModel, Outcome, dataframe, 
categorizationType=c("Raw","Categorical","ZCategorical",
"RawZCategorical","RawZTail","RawTail"), 
FitType = c("LOGIT", "LM", "COX"), 
rankingTest = c("zIDI", "zNRI","IDI", "NRI","NeRI","Ztest","AUC","CStat","Kendall"), 
cateGroups = c(0.1, 0.9),
raw.dataFrame=NULL,
description=".",
uniType=c("Binary","Regression")
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varList}{
  A data frame with the candidate variables to be ranked
%%     ~~Describe \code{varList} here~~
}
  \item{baseModel}{
  A base model containing the covariates ie: Gender, BMI, Height to be adjusted
%%     ~~Describe \code{baseModel} here~~
}
  \item{Outcome}{
  The binary outcome 

}
  \item{dataframe}{
  The data frame with column features and row observations

}
  \item{categorizationType}{
  How variables will be ranked: 
  Raw: As given in the data frame. 
  Categorical: variables are broken into the p.values categories.
  ZCategorical: Categories will be weighted by it z value.
  RawZCategorical: Categories plus raw.
  RawZTail: The raw variables plus the two tails.
%%     ~~Describe \code{categorizationType} here~~
}
  \item{FitType}{
  The fitting function to be used:
  LM: Lineal regression fitted by gls.
  LOGIT: Logistic regression fitted by gls.
  COX: cox proportional model fitted by coxph (Survival package).
%%     ~~Describe \code{FitType} here~~
}
  \item{rankingTest}{
  The test to be used to rank the variables.
  IDI: Integrated Discriminant Index.
  NRI: Net reclassification Improvement.
  NeRI: Net Residual Improvement.
  Z: Z value of the model fit.
  AUC: ROC AUC
%%     ~~Describe \code{rankingTest} here~~
}
  \item{cateGroups}{
  A list of p.values that will be used to create the categories.
%%     ~~Describe \code{cateGroups} here~~
}
  \item{raw.dataFrame}{
 The dataframe with the original (non-transformed) data
}
  \item{description}{
  The varList's column name that contains the variable description
}
  \item{uniType}{
  The type of univariate analysis: Binary classification or Regression
}

}
\details{
This univariate function will create valid dummy categorical variables if and only if, the dataframe has been z-standarized.
The provided p.values for categorization will be converted to its corresponding z-value and this z-value will be used for the 
variable categorization. If the data frame has not been z-standardized the categorization analysis provided by the function is wrong.

}
\value{
Returns a sorted frame of the univariate analysis.
For Binary Classification it will return a data frame with the following columns: \cr 
Name: The raw variable or dummy variable if the data has been categorized \cr
parent: The raw variable name \cr
descrip: If the description is present, it will describe the parent variable\cr
cohortMean: The entire set mean\cr
cohortStd: the entire set standard deviation\cr
cohortKSD: the entire set D value of the Kolmogorov-Smirnov(KS) test compared to a normal distribution\cr
cohortKSP: The associated p-value of the test compared\cr
caseMean: The mean value of the cases (Outcome == 1)\cr
caseStd: The standard deviation value of the cases (Outcome == 1)\cr
caseKSD: The (KS) D-statistics of the cases (Outcome == 1)\cr
caseKSP: The (KS) p-value of the cases (Outcome == 1)\cr
caseZKSD: The (KS) D-statistics of the z-standarized cases (Outcome == 1)\cr
caseZKSP: The (KS) p-value of the z-standarized cases (Outcome == 1)\cr
controlMean: he mean value of the cases (Outcome == 0)\cr
controlStd: The standard deviation value of the cases (Outcome == 0)\cr
controlKSD: The (KS) D-statistics of the cases (Outcome == 0)\cr
controlKSP: The (KS) p-value of the cases (Outcome == 0)\cr
controlZKSD: The (KS) D-statistics of the z-standarized cases (Outcome == 0)\cr
controlZKSP: The (KS) p-value of the z-standarized cases (Outcome == 0)\cr
t.Rawvalue: The normal inverse p-value (Z-value) of the t-test performed on the raw data set\cr
t.Zvalue: The normal inverse p-value (Z-value) of the t-test performed on the dataframe\cr
wilcox.Zvalue: The normal inverse p-value (Z-value) of the Mann-Whitney test performed on the dataframe\cr
ZGLM:: The z-value returned by the glm or COX procedure of the tested adjusted variable\cr
zNRI: The z-value returned by the improveProb procedure (HMISC) for the NRI\cr
zIDI: The z-value returned by the improveProb procedure (HMISC) for the IDI\cr
zNeRI: The z-value returned by the improvedResiduals procedure (FRESA.CAD) for the NeRI\cr
ROCAUC: The area under the curve returned by the roc procedure (pROC) \cr
cStatCorr: The c-stat correlation coefficient returned by the rcorr.cens procedure (HMISC)\cr
NRI: The NRI statistics: improveProb procedure (HMISC) \cr
IDI: The IDI statistics: improveProb procedure (HMISC) \cr
NeRI: The NeRI statistics: improvedResiduals procedure (FRESA.CAD)\cr
kendall.r: The Kendall tau correlation to the binary outcome\cr
kendall.p: The p-value of the Kendall tau correlation to the binary outcome\cr
TstudentRes.p: The p-value of the improvement in residuals as evaluated by the paired t-test\cr
WilcoxRes.p: The p-value of the improvement in residuals as evaluated by the paired Mann-Whitney test\cr
FRes.p: The p-value of the improvement in residuals variance as evaluated by the F test \cr
caseN_Z_Low_Tail: The number of case subjects in the low tail\cr
caseN_Z_Hi_Tail: The number of case subjects in the top tail\cr
controlN_Z_Low_Tail: The number of control subjects in the low tail\cr
controlN_Z_Hi_Tail: The number of control subjects in the top tail\cr

For residual analysis it will return a data frame with the following columns
Name: The raw variable or dummy variable if the data has been categorized\cr
parent: The raw variable name\cr
descrip: If the description is present, it will describe the parent variable\cr
cohortMean: The entire set mean\cr
cohortStd: the entire set standard deviation\cr
cohortKSD: the entire set D value of the Kolmogorov-Smirnov(KS) test compared to a normal distribution\cr
cohortKSP: The associated p-value of the test compared to a normal distribution\cr
cohortZKSD: The entire set D value of the Kolmogorov-Smirnov(KS) test of the z standardized data compared to a normal distribution\cr
cohortZKSP: The associated p-value of the KS test of the z standardized data compared to a normal distribution\cr
ZGLM:: The z-value returned by the glm or COX procedure of the tested adjusted variable\cr
zNeRI: The z-value returned by the improvedResiduals procedure (FRESA.CAD) for the NeRI\cr
NeRI: The NeRI statistics: improvedResiduals procedure (FRESA.CAD)\cr
cStatCorr: The c-stat correlation coefficient returned by the rcorr.cens procedure (HMISC)\cr
spearman.r: The Spearman rho correlation between the outcome and the variable\cr
pearson.r: The Pearson rho correlation between the outcome and the variable\cr
kendall.r: The Kendaall Tau correlation between the outcome and the variable\cr
kendall.p: The Kendaall p.value of the Tau correlation between the outcome and the variable\cr
TstudentRes.p: The p-value of the improvement in residuals as evaluated by the paired t-test\cr
WilcoxRes.p: The p-value of the improvement in residuals as evaluated by the paired Mann-Whitney test\cr
FRes.p: The p-value of the improvement in residuals variance as evaluated by the F test \cr
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
%%% Add one or more standard keywords, see file 'KEYWORDS' in the
%%% R documentation directory.
%%}
\keyword{Data_Inspection}

