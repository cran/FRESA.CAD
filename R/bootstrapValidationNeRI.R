bootstrapValidationNeRI <-
function(fraction=1.00,loops=200,model.formula,Outcome,data,type=c("LM","LOGIT","COX"),plots=TRUE)
{



    type <- match.arg(type);
	basemodel <- NULL
	if (class(model.formula) != "formula")
	{
		model.formula = formula(model.formula);
	}
	varsList <- as.list(attr(terms(model.formula),"variables"))
	if (length(varsList)>2)
	{
		outn <- length(table(data[,Outcome]))
		basemodel <- modelFitting(model.formula,data,type)
		
		if ( inherits(basemodel, "try-error"))
		{
			stop("BootstapValidation: Initial model Error\n")
		}

		modelMat <- model.matrix(model.formula,data);
		if (type=="COX")
		{
			response <- data[,all.vars(model.formula)[1:2]];
		}
		else
		{
			response <- data[,Outcome];
		}

		 output<-.Call("bootstrapValidationNeRICpp", fraction, loops, modelMat,type,data.matrix(response),1);


		bootmodel <- basemodel;
		bootmodel$coefficients <-  apply(output$bcoef,2,mean,trim = 0.25, na.rm = TRUE)
		names(bootmodel$coefficients) <- names(basemodel$coefficients);

		pr <- predictForFresa(bootmodel,testData=data,predictType = 'linear');
		bootmodel$linear.predictors <- pr;
		colnames(output$NeRi) <- attr(terms(model.formula),'term.labels');
		colnames(output$tpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$wpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$spvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$Fpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$test_tpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$test_wpvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$test_spvalue) <- attr(terms(model.formula),'term.labels');
		colnames(output$test_Fpvalue) <- attr(terms(model.formula),'term.labels');
		testRMSE <- sqrt(mean(output$testResiduals^2));
		trainRMSE <- sqrt(mean(output$trainResiduals^2));
		

		if (plots)
		{
			if (outn>2) 
			{	
				par(mfrow=c(1,3))
				plot(output$testPrediction ~ output$testOutcome,main="Prediction~Outcome");
				plot(output$testResiduals ~ output$testOutcome,main="Residuals");
			}
			else
			{
				par(mfrow=c(2,2))
				pROC::roc(output$testOutcome, output$testPrediction,auc=TRUE,plot=TRUE,smooth=FALSE,print.auc=TRUE);
				boxplot(output$testPrediction ~ output$testOutcome,main="Prediction~Outcome");
				boxplot(output$testResiduals ~ output$testOutcome,main="Residuals");
			}
			plot(ecdf(output$trainSampledRMSE),main="RMSE",col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
			abline(v=testRMSE,col = "red");
			par(mfrow=c(1,1))
			cat("Train RMSE :",trainRMSE," Test RMSE :",testRMSE," F Ratio :",(testRMSE/trainRMSE)^2,"Original RMSE :",output$startRMSE," Boot RMSE :",output$bootRMSE," F Ratio :",(output$bootRMSE/output$startRMSE)^2,"\n");
		}

		result <- structure(llist( data = data,
			outcome=data[,Outcome],
			boot.model=bootmodel,
			NeRis = output$NeRi,
			tStudent.pvalues = output$tpvalue,
			wilcox.pvalues = output$wpvalue,
			bin.pvlaues = output$spvalue,
			F.pvlaues = output$Fpvalue,
			test.tStudent.pvalues = output$test_tpvalue,
			test.wilcox.pvalues = output$test_wpvalue,
			test.bin.pvlaues = output$test_spvalue,
			test.F.pvlaues = output$test_Fpvalue,
			testPrediction=output$testPrediction,
			testOutcome=output$testOutcome,
			testResiduals=output$testResiduals,
			trainPrediction=output$trainPrediction,
			trainOutcome=output$trainOutcome,
			trainResiduals=output$trainResiduals,
			testRMSE=testRMSE,
			trainRMSE=trainRMSE,
			trainSampledRMSE=output$trainSampledRMSE,
			testSampledRMSE=output$testSampledRMSE
		), class = "bootstrapValidationNeRI");


		return (result);
	}
	else
	{
		return (NULL);
	}
}
