getVarNeRI <-
function (object,data,Outcome="Class", type=c("LM","LOGIT","COX"),testData=NULL) 
{
	if (is.null(testData))
	{
		testData <- data;
	}
  
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
	fullModel <- modelFitting(ftmp,data,type);

	fullResiduals <- residualForNeRIs(fullModel,data,Outcome);
	testResiduals <- residualForNeRIs(fullModel,testData,Outcome);

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
			redModel <- modelFitting(ftmp,data,type)

			if ( inherits(redModel, "try-error"))
			{
				redModel <- fullModel;
			}
			
			redResiduals <- residualForNeRIs(redModel,data,Outcome);
			redTestResiduals <- residualForNeRIs(redModel,testData,Outcome);
			iprob <- improvedResiduals(redResiduals,fullResiduals);
			testiprob <- improvedResiduals(redTestResiduals,testResiduals);
			model_tpvalue <- append(model_tpvalue,iprob$tP.value);
			model_bpvalue <- append(model_bpvalue,iprob$BinP.value);
			model_wpvalue <- append(model_wpvalue,iprob$WilcoxP.value);
			model_fpvalue <- append(model_fpvalue,iprob$FP.value);
			model_neri <- append(model_neri,iprob$NeRI);
			testmodel_tpvalue <- append(testmodel_tpvalue,testiprob$tP.value);
			testmodel_bpvalue <- append(testmodel_bpvalue,testiprob$BinP.value);
			testmodel_wpvalue <- append(testmodel_wpvalue,testiprob$WilcoxP.value);
			testmodel_fpvalue <- append(testmodel_fpvalue,testiprob$FP.value);
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
