medianPredict <-
function (formulaList,traindata,newdata=NULL, predictType = c("prob", "linear"),modelType = c("LOGIT", "LM","COX"),Outcome="CLASS",k=0,...) 
{

	if (k >= 0)
	{
		casesample = subset(traindata,get(Outcome)  == 1);
		controlsample = subset(traindata,get(Outcome) == 0);
		
		minTrainSamples <- min(nrow(casesample),nrow(controlsample));
		
		KnnTrainSet <- rbind(casesample[sample(1:nrow(casesample),minTrainSamples,replace=FALSE),],controlsample[sample(1:nrow(controlsample),minTrainSamples,replace=FALSE),])
	}


	if (k==0)
	{
		k = 2*as.integer(sqrt(minTrainSamples/2)) + 1;
	}
	if (is.null(newdata))
	{
		newdata <- traindata
	}

	ftmp <- formula(formulaList[1])
	bestmodel <- modelFitting(ftmp,traindata,modelType)	
	out <- predictForFresa(bestmodel,newdata,predictType);
	if (k>0)
	{
		outKNN <- getKNNpredictionFromFormula(ftmp,KnnTrainSet,newdata,outcome=Outcome,k)$binProb
	}
	else
	{
		outKNN <- NULL;
		medianKNN <- NULL;
	}
	
	for (i in 2:length(formulaList))
	{
		ftmp <- formula(formulaList[i]);
		out <- cbind(out,predictForFresa(modelFitting(ftmp,traindata,modelType),newdata,predictType,...));
		if (k>0) 
		{
			outKNN <- cbind(outKNN,getKNNpredictionFromFormula(ftmp,KnnTrainSet,newdata,outcome=Outcome,k)$binProb);
		}
	}
	medianout <- rowMedians(out);
	if (k>0) 
	{
		medianKNN <- rowMedians(outKNN);
	}
	result <- list(modelPredict=medianout,
	KNN.Predict=medianKNN,predictions=out,KNNpredictions=outKNN)
    return (result)
}


