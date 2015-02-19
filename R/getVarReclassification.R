getVarReclassification <-
function (object,data,Outcome="Class", type = c("LOGIT", "LM","COX"),testData=NULL) 
{

	if (is.null(testData))
	{
		testData <- data;
	}

	model_zidi <- vector();
	model_znri <- vector();
	model_idi <- vector();
	model_nri <- vector();

	t.model_zidi <- vector();
	t.model_znri <- vector();
	t.model_idi <- vector();
	t.model_nri <- vector();
  
	varsList <- as.list(attr(terms(object),"variables"))

		
	outCome = paste(varsList[2]," ~ 1");
	startlist = 3;
	frm1 = outCome;
	for ( i in startlist:length(varsList))
	{
		frm1 <- paste(frm1,paste(" + ",varsList[i]));
	}
	ftmp <- formula(frm1);
	
	fullModel <- modelFitting(ftmp,data,type)
	if ( !inherits(fullModel, "try-error"))
	{
	
		fullPredict_train <- predictForFresa(fullModel,data,'prob');
		fullPredict <- predictForFresa(fullModel,testData, 'prob');

		for ( i in startlist:length(varsList))
		{
		
			frm1 = outCome;
				for ( j in startlist:length(varsList))
				{
					if (i!=j)
					{
						frm1 <- paste(frm1,paste(" + ",varsList[j]));
					}
				}
			ftmp <- formula(frm1);
			redModel <- modelFitting(ftmp,data,type)
			if (inherits(redModel, "try-error"))
			{
				redModel <- fullModel;
			}

			redPredict <- predictForFresa(redModel,testData,'prob');
			redPredict_train <- predictForFresa(redModel,data,'prob');
			
			iprob <- improveProb(redPredict,fullPredict,testData[,Outcome]);
			iprob_t <- improveProb(redPredict_train,fullPredict_train,data[,Outcome]);


			model_zidi <- append(model_zidi,iprob$z.idi);
			model_idi <- append(model_idi,iprob$idi);
			model_nri <- append(model_nri,iprob$nri);
			model_znri <- append(model_znri,iprob$z.nri);

			t.model_zidi <- append(t.model_zidi,iprob_t$z.idi);
			t.model_idi <- append(t.model_idi,iprob_t$idi);
			t.model_nri <- append(t.model_nri,iprob_t$nri);
			t.model_znri <- append(t.model_znri,iprob_t$z.nri);
		}

	}
	
	 result <- list(
	 z.IDIs=t.model_zidi,
	 z.NRIs=t.model_znri,
	 IDIs=t.model_idi,
	 NRIs=t.model_nri,
	 testData.z.IDIs=model_zidi,
	 testData.z.NRIs=model_znri,
	 testData.IDIs=model_idi,
	 testData.NRIs=model_nri
	 );

    return (result)
}
