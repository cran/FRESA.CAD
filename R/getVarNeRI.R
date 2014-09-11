getVarNeRI <-
function (object,dataframe,Outcome="Class", type=c("LM","LOGIT","COX"),testdata=NULL) 
{
	if (is.null(testdata))
	{
		testdata <- dataframe;
	}
    type <- match.arg(type);
  
	varsList <- as.list(attr(terms(object),"variables"))
	termList <- as.list(attr(terms(object),"term.labels"))
	
	outCome = paste(varsList[2]," ~ 1");
	startlist = 3;
	frm1 = outCome;
	for ( i in startlist:length(varsList))
	{
		frm1 <- paste(frm1,paste(" + ",varsList[i]));
	}
	ftmp <- formula(frm1);
	fullModel <- modelFitting(ftmp,dataframe,type);

	fullResiduals <- residualForNeRIs(fullModel,newdata=dataframe,Outcome);
	testResiduals <- residualForNeRIs(fullModel,newdata=testdata,Outcome);

	model_tpvalue <- vector();
	model_bpvalue <- vector();
	model_neri <- vector();
	model_wpvalue <- vector();
	model_fpvalue <- vector();
	testmodel_tpvalue <- vector();
	testmodel_bpvalue <- vector();
	testmodel_neri <- vector();
	testmodel_wpvalue <- vector();
	testmodel_fpvalue <- vector();
	if (length(termList)>0)
	{
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

			if ( inherits(redModel, "try-error"))
			{
				redModel <- fullModel;
			}
			
			redResiduals <- residualForNeRIs(redModel,newdata=dataframe,Outcome);
			redTestResiduals <- residualForNeRIs(redModel,newdata=testdata,Outcome);
			iprob <- improvedResiduals(redResiduals,fullResiduals);
			testiprob <- improvedResiduals(redTestResiduals,testResiduals);
			model_tpvalue <- append(model_tpvalue,iprob$t.test.pValue);
			model_bpvalue <- append(model_bpvalue,iprob$binom.pValue);
			model_wpvalue <- append(model_wpvalue,iprob$wilcox.pValue);
			model_fpvalue <- append(model_fpvalue,iprob$F.test.pValue);
			model_neri <- append(model_neri,iprob$NeRI);
			testmodel_tpvalue <- append(testmodel_tpvalue,testiprob$t.test.pValue);
			testmodel_bpvalue <- append(testmodel_bpvalue,testiprob$binom.pValue);
			testmodel_wpvalue <- append(testmodel_wpvalue,testiprob$wilcox.pValue);
			testmodel_fpvalue <- append(testmodel_fpvalue,testiprob$F.test.pValue);
			testmodel_neri <- append(testmodel_neri,testiprob$NeRI);
		}
	}
	
	 result <- list(
	 tP.value=model_tpvalue,
	 BinP.value=model_bpvalue,
	 WilcoxP.value=model_wpvalue,
	 FP.value=model_fpvalue,
	 NeRIs=model_neri,
	 testData.tP.value=testmodel_tpvalue,
	 testData.BinP.value=testmodel_bpvalue,
	 testData.WilcoxP.value=testmodel_wpvalue,
	 testData.FP.value=testmodel_fpvalue,
	 testData.NeRIs=testmodel_neri
	 );

    return (result)
}
