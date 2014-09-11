getKNNpredictionFromFormula <-
function (modelFormula,trainData,testData,outcome="CLASS",k=3) 
{
if (!require(class)) {
  knn <- function(...) {}
  install.packages("class", dependencies = TRUE)
  library(class)
}


	temslist <- attr(terms(modelFormula),"term.labels")
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

	# for (i in 1:length(temslist))
	# {
		# varlist <- append(varlist,str_replace_all(unlist(strsplit(str_replace_all(str_replace_all(temslist[i],"I\\("," "),"\\("," "),"\\*"))[1]," ",""))
	# }
	# varlist <- as.vector(rownames(table(varlist)))

	print(varlist)
		
	knnclass <- knn(as.data.frame(trainData[,varlist]),as.data.frame(testData[,varlist]),factor(trainData[,outcome]),k,prob=TRUE)
	prop <- attributes(knnclass)
	binProp <-abs(prop$prob-1*(knnclass=="0"))
	result <- list(prediction=knnclass,
	prob=prop,
	binProb=binProp,
	featureList=varlist)
    return (result)
}
