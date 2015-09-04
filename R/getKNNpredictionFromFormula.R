getKNNpredictionFromFormula <-
function (model.formula,trainData,testData,Outcome="CLASS",nk=3) 
{

if (!requireNamespace("class", quietly = TRUE)) {
   install.packages("class", dependencies = TRUE)
} 


	temslist <- attr(terms(model.formula),"term.labels")
	varlist <- vector()

	for (i in 1:length(temslist))
	{
		varlist <- append(varlist,str_replace_all(unlist(strsplit(
						str_replace_all(
							str_replace_all(
								str_replace_all(
									str_replace_all(temslist[i],"I\\("," ")
								,"\\("," ")
							,">","\\*")
						,"<","\\*")
				,"\\*"))[1]," ",""))
	}
	varlist <- as.vector(rownames(table(varlist)))

		
	knnclass <- class::knn(as.data.frame(trainData[,varlist]),as.data.frame(testData[,varlist]),factor(trainData[,Outcome]),nk,prob=TRUE)
	prop <- attributes(knnclass)
	binProp <-abs(prop$prob-1*(knnclass=="0"))
	result <- list(prediction=knnclass,
	prob=prop,
	binProb=binProp,
	featureList=varlist)
    return (result)
}
