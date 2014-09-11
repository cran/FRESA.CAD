getVarReclassification <-
function (object,dataframe,Outcome="Class", type = c("LOGIT", "LM","COX"),independentFrame=NULL) 
{
    type <- match.arg(type);

	if (is.null(independentFrame))
	{
		independentFrame <- dataframe;
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
	
	fullModel <- modelFitting(ftmp,dataframe,type)
	if ( !inherits(fullModel, "try-error"))
	{
	
		fullPredict_train <- predictForFresa(fullModel,newdata=dataframe,type = 'prob');
		fullPredict <- predictForFresa(fullModel,newdata=independentFrame,type = 'prob');

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
			redModel <- modelFitting(ftmp,dataframe,type)
			if (inherits(redModel, "try-error"))
			{
				redModel <- fullModel;
			}

			redPredict <- predictForFresa(redModel,newdata=independentFrame,type = 'prob');
			redPredict_train <- predictForFresa(redModel,newdata=dataframe,type = 'prob');
			
			iprob <- improveProb(redPredict,fullPredict,independentFrame[,Outcome]);
			iprob_t <- improveProb(redPredict_train,fullPredict_train,dataframe[,Outcome]);


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
	 z.IDIs=model_zidi,
	 z.NRIs=model_znri,
	 IDIs=model_idi,
	 NRIs=model_nri,
	 tz.IDIs=t.model_zidi,
	 tz.NRIs=t.model_znri,
	 tIDIs=t.model_idi,
	 tNRIs=t.model_nri
	 );

    return (result)
}
