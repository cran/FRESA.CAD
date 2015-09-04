medianPredict <-
function (formulaList,trainData,testData=NULL, predictType = c("prob", "linear"),type = c("LOGIT", "LM","COX"),Outcome="CLASS",nk=0,...) 
{

	if (nk >= 0)
	{
		casesample = subset(trainData,get(Outcome)  == 1);
		controlsample = subset(trainData,get(Outcome) == 0);
		
		minTrainSamples <- min(nrow(casesample),nrow(controlsample));
		
		KnnTrainSet <- rbind(casesample[sample(1:nrow(casesample),minTrainSamples,replace=FALSE),],controlsample[sample(1:nrow(controlsample),minTrainSamples,replace=FALSE),])
	}


	if (nk==0)
	{
		nk = 2*as.integer(sqrt(minTrainSamples/2)) + 1;
	}
	if (is.null(testData))
	{
		testData <- trainData
	}

	ftmp <- formula(paste(Outcome,"~",formulaList[[1]]))
	bestmodel <- modelFitting(ftmp,trainData,type)	
	out <- predictForFresa(bestmodel,testData,predictType);
	if (nk>=0)
	{
		outKNN <- getKNNpredictionFromFormula(ftmp,KnnTrainSet,testData,Outcome=Outcome,nk)$binProb
	}
	else
	{
		outKNN <- NULL;
		medianKNN <- NULL;
	}
	
	for (i in 2:length(formulaList))
	{
		ftmp <- formula(paste(Outcome,"~",formulaList[[i]]));
		out <- cbind(out,predictForFresa(modelFitting(ftmp,trainData,type),testData,predictType,...));
		if (nk>=0) 
		{
			outKNN <- cbind(outKNN,getKNNpredictionFromFormula(ftmp,KnnTrainSet,testData,Outcome=Outcome,nk)$binProb);
		}
	}
	medianout <- rowMedians(out);
	if (nk>=0) 
	{
		medianKNN <- rowMedians(outKNN);
	}
	result <- list(medianPredict=medianout,
	medianKNNPredict=medianKNN,predictions=out,KNNpredictions=outKNN)
    return (result)
}