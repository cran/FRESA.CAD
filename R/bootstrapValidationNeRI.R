bootstrapValidationNeRI <-
function(fraction=1.00,loops=200,model.formula,Outcome,dataframe,type=c("LM","LOGIT","COX"),plots=TRUE)
{
    type <- match.arg(type);
  
	datasize = nrow(dataframe);
	totSamples = as.integer(fraction*datasize);
	basemodel <- NULL
	 
#	cat("In NeRI Bootstrap \n")
	basemodel <- modelFitting(model.formula,dataframe,type)

	startRSME <- sqrt(mean(residualForNeRIs(basemodel,dataframe,Outcome)^2))

	if ( inherits(basemodel, "try-error"))
	{
		stop("BootstapValidation: Initial model Error\n")
	}
#	print(summary(basemodel))

	testOutcome <- vector();
	testPrediction <- vector();
	testResiduals <- vector();
	trainOutcome <- vector();
	trainPrediction <- vector();
	trainResiduals <- vector();
	testSampledRSME <- vector();
	trainSampleRSME <- vector();
	
	outn <- length(table(dataframe[,Outcome]))
	valinserted = 0;
	
	for ( doOver in 1:loops)
	{ 
		trainSamples <- sample(1:datasize, totSamples, replace=TRUE);
		mysampleCases <- dataframe[trainSamples,];
		testCases <-dataframe[-trainSamples,];
		testResampledCases <-testCases[sample(1:nrow(testCases), totSamples, replace=TRUE),];
		
		model <- modelFitting(model.formula,mysampleCases,type)
		
		if ( !inherits(model, "try-error"))
		{
			cf <- as.vector(model$coefficients);
			
			residual <- residualForNeRIs(model,testCases,Outcome);
			testOutcome <- append(testOutcome,testCases[,Outcome]);
			testPrediction <- append(testPrediction,testCases[,Outcome]+residual);
			testResiduals <- append(testResiduals,residual);
			testSampledRSME <- append(testSampledRSME,sqrt(mean(residual^2)));

			residual <- residualForNeRIs(model,mysampleCases,Outcome);
			trainOutcome <- append(trainOutcome,mysampleCases[,Outcome]);
			trainPrediction <- append(trainPrediction,mysampleCases[,Outcome]+residual);
			trainResiduals <- append(trainResiduals,residual);
			trainSampleRSME <- append(trainSampleRSME,sqrt(mean(residual^2)));

			modelReclas <- getVarNeRI(model,dataframe=mysampleCases,Outcome=Outcome,type=type,testdata=testResampledCases);
			valinserted = valinserted + 1;

			if (valinserted == 1) 
			{
				bcoef <- t(matrix(cf));
				NeRi <- t(matrix(modelReclas$NeRIs));
				tpvalue <- t(matrix(modelReclas$tP.value));
				wpvalue <- t(matrix(modelReclas$WilcoxP.value));
				spvalue <- t(matrix(modelReclas$BinP.value));
				Fpvalue <- t(matrix(modelReclas$FP.value));
				test.tpvalue <- t(matrix(modelReclas$testData.tP.value));
				test.wpvalue <- t(matrix(modelReclas$testData.WilcoxP.value));
				test.spvalue <- t(matrix(modelReclas$testData.BinP.value));
				test.Fpvalue <- t(matrix(modelReclas$testData.FP.value));
#				print(modelReclas$BinP.value)
#				print(modelReclas$testData.BinP.value)
#				print(spvalue)
			}
			else 
			{
				bcoef <- rbind(bcoef,t(matrix(cf)));	  
				NeRi <- rbind(NeRi,t(matrix(modelReclas$NeRIs)));
				tpvalue <- rbind(tpvalue,t(matrix(modelReclas$tP.value)));
				wpvalue <- rbind(wpvalue,t(matrix(modelReclas$WilcoxP.value)));
				spvalue <- rbind(spvalue,t(matrix(modelReclas$BinP.value)));
				Fpvalue <- rbind(Fpvalue,t(matrix(modelReclas$FP.value)));
				test.tpvalue <- rbind(test.tpvalue,t(matrix(modelReclas$testData.tP.value)));
				test.wpvalue <- rbind(test.wpvalue,t(matrix(modelReclas$testData.WilcoxP.value)));
				test.spvalue <- rbind(test.spvalue,t(matrix(modelReclas$testData.BinP.value)));
				test.Fpvalue <- rbind(test.Fpvalue,t(matrix(modelReclas$testData.FP.value)));

			}
			if (plots && ((doOver %% 20) == 0))
			{
				testRSME <- sqrt(mean(testResiduals^2));
				trainRSME <- sqrt(mean(trainResiduals^2));
				cat("Train RMSE :",trainRSME," Test RMSE :",testRSME," F Ratio :",(testRSME/trainRSME)^2,"\n");
				if (outn>2) 
				{	
					par(mfrow=c(1,3))
					plot(testPrediction ~ testOutcome,main="Prediction~Outcome");
					plot(testResiduals ~ testOutcome,main="Residuals");
				}
				else
				{
					par(mfrow=c(2,2))
					roc(testOutcome, testPrediction,auc=TRUE,plot=TRUE,smooth=FALSE,print.auc=TRUE);
					boxplot(testPrediction ~ testOutcome,main="Prediction~Outcome");
					boxplot(testResiduals ~ testOutcome,main="Residuals");
				}
				plot(ecdf(trainSampleRSME),main="RMSE",col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
				abline(v=testRSME,col = "red");
			}
		}
		else
		{
			cat(doOver,": Warning Bootstrap fitting error \n");
#			print(trainmodel$coefficients);
		}


	}

	bootmodel <- basemodel;
	bootmodel$coefficients <- colMedians(bcoef,na.rm = TRUE);
	
	names(bootmodel$coefficients) <- names(basemodel$coefficients);
	pr <- predictForFresa(bootmodel,newdata=dataframe,type = 'linear');
   	bootmodel$linear.predictors <- pr;


#	print(summary(bootmodel));
	

    colnames(NeRi) <- attr(terms(model.formula),'term.labels');
	colnames(tpvalue) <- attr(terms(model.formula),'term.labels');
	colnames(wpvalue) <- attr(terms(model.formula),'term.labels');
	colnames(spvalue) <- attr(terms(model.formula),'term.labels');
	colnames(Fpvalue) <- attr(terms(model.formula),'term.labels');
	colnames(test.tpvalue) <- attr(terms(model.formula),'term.labels');
	colnames(test.wpvalue) <- attr(terms(model.formula),'term.labels');
	colnames(test.spvalue) <- attr(terms(model.formula),'term.labels');
	colnames(test.Fpvalue) <- attr(terms(model.formula),'term.labels');


	colnames(bcoef) <- names(basemodel$coefficients);
	

	testRSME <- sqrt(mean(testResiduals^2));
	trainRSME <- sqrt(mean(trainResiduals^2));
	bootRSME <- sqrt(mean(residualForNeRIs(bootmodel,dataframe,Outcome)^2))
	cat("Train RMSE :",trainRSME," Test RMSE :",testRSME," F Ratio :",(testRSME/trainRSME)^2,"Original RMSE :",startRSME," Boot RMSE :",bootRSME," F Ratio :",(bootRSME/startRSME)^2,"\n");

	if (outn>2) 
	{	
		par(mfrow=c(1,3))
		plot(testPrediction ~ testOutcome,main="Prediction~Outcome");
		plot(testResiduals ~ testOutcome,main="Residuals");
	}
	else
	{
		par(mfrow=c(2,2))
		roc(testOutcome, testPrediction,auc=TRUE,plot=TRUE,smooth=FALSE,print.auc=TRUE);
		boxplot(testPrediction ~ testOutcome,main="Prediction~Outcome");
		boxplot(testResiduals ~ testOutcome,main="Residuals");
	}
	plot(ecdf(trainSampleRSME),main="RMSE",col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
	abline(v=testRSME,col = "red");
	par(mfrow=c(1,1))

	result <- structure(llist(outcome=dataframe[,Outcome],
		boot.model=bootmodel,
		NeRis = NeRi,
		tStudent.pvalues = tpvalue,
		wilcox.pvalues = wpvalue,
		bin.pvlaues = spvalue,
		F.pvlaues = Fpvalue,
		test.tStudent.pvalues = test.tpvalue,
		test.wilcox.pvalues = test.wpvalue,
		test.bin.pvlaues = test.spvalue,
		test.F.pvlaues = test.Fpvalue,
		testPrediction=testPrediction,
		testOutcome=testOutcome,
		testResiduals=testResiduals,
		trainPrediction=testPrediction,
		trainOutcome=testOutcome,
		trainResiduals=testResiduals,
		testRMSE=testRSME,
		trainRMSE=trainRSME,
		trainSampleRSME=trainSampleRSME,
		testSampledRSME=testSampledRSME
	), class = "bootstrapValidationNeRI");
	
	return (result);
}
