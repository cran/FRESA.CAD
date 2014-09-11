bootstrapValidation <-
function(fraction=1.00,loops=200,model.formula,Outcome,dataframe,type=c("LM","LOGIT","COX"),plots=TRUE)
{
      type <- match.arg(type);
  
	 casesample = subset(dataframe,get(Outcome)  == 1);
	 controlsample = subset(dataframe,get(Outcome) == 0);

	 sizecases = nrow(casesample);
	 sizecontrol = nrow(controlsample);
	 minsize = min(sizecases,sizecontrol);
	 
	 totSamples = as.integer(fraction*minsize);
	 
	 accuracy <- vector();
	 sensitivity <- vector();
	 specificity <- vector();

	 taccuracy <- vector();
	 tsensitivity <- vector();
	 tspecificity <- vector();
	 trainRoc <- vector();
	 
	basemodel <- modelFitting(model.formula,dataframe,type)


     basepredict <- predictForFresa(basemodel,newdata=dataframe,type = 'linear');
	 basemodel$linear.predictors <- basepredict;

	 framesize = nrow(dataframe);
	 acc = 0.0;
	 sen = 0.0;
	 spe = 0.0;
	 for (i in 1:framesize)
	 {
		if ((dataframe [i,Outcome] > 0) && (basepredict[i] > 0) ) 
		{
			acc = acc + 1.0; 
			sen = sen + 1.0;
		}
		if ((dataframe [i,Outcome] == 0) && (basepredict[i] < 0) ) 
		{
			acc = acc + 1.0; 
			spe = spe + 1.0;
		}
	}
	baseAcc = acc/framesize;
	baseSen = sen/sizecases;
	baseSpe = spe/sizecontrol;
	gacc = 0; 
	gsen = 0; 
	gspe = 0; 
	smaptot = 0;
	smapsen = 0;
	smapspe = 0;
	testoutcome <- vector();
	testprediction <- vector();
	
	sumwtdcf <- rep(0,length(basemodel$coefficients));
	sumwts <- sumwtdcf;
	valinserted = 0

	for ( doOver in 1:loops)
	{ 
		rsamp = totSamples;
		caseSamp <- sample(1:sizecases, rsamp, replace=TRUE);
		contSamp <- sample(1:sizecontrol, rsamp, replace=TRUE);
		mysampleCases = casesample[caseSamp,];
		mysampleControl = controlsample [contSamp,];
		  
		myTestCases = casesample[-caseSamp,];
		myTestControl = controlsample [-contSamp,];
		
#		leftCases = rsamp - nrow(myTestCases)
#		leftControl = rsamp - nrow(myTestControl)

		trainingSample = rbind(mysampleCases,mysampleControl);
		testSample = rbind(myTestCases,myTestControl);
		
		indsampleCases <- myTestCases[sample(1:nrow(myTestCases), rsamp, replace=TRUE),];
		indsampleControl <- myTestControl[sample(1:nrow(myTestControl), rsamp, replace=TRUE),];

		indpendentSample = rbind(indsampleCases,indsampleControl);

		trainmodel <- modelFitting(model.formula,trainingSample,type)
		
		if ( !inherits(trainmodel, "try-error"))
		{
		  

			  trainmodel$linear.predictors <- predictForFresa(trainmodel,newdata=trainingSample,type = 'linear');
			  p <- predictForFresa(trainmodel,newdata=testSample,type = 'linear');

			  aucTrain <- roc(trainingSample[,Outcome], trainmodel$linear.predictors,auc=TRUE,plot=FALSE,smooth=FALSE)$auc

			  trainRoc <- append(trainRoc,aucTrain);

			  testsize = nrow(testSample);
			  testsizeCase = nrow(myTestCases);
			  testsizeCont = nrow(myTestControl);

			  sen = sum( 1*((testSample [,Outcome] > 0)*( p >= 0 )) , na.rm = TRUE);
			  spe = sum( 1*((testSample [,Outcome] == 0)*( p < 0 )) , na.rm = TRUE);
			  acc = sen+spe;
		  
			  plog <- predictForFresa(trainmodel,newdata=testSample,type = 'prob');
			  testoutcome <- append(testoutcome,testSample [,Outcome]);
			  testprediction <- append(testprediction,plog);
			  
			  gacc = gacc + acc;
			  smaptot = smaptot+testsize;
			  gsen = gsen + sen; 
			  smapsen = smapsen + testsizeCase;
			  gspe = gspe + spe; 
			  smapspe = smapspe + testsizeCont
			  
			  acc = acc/testsize;
			  sen = sen/testsizeCase;
			  spe = spe/testsizeCont;

			  
			  accuracy <- append(accuracy,acc);
			  sensitivity <- append(sensitivity,sen);
			  specificity <- append(specificity,spe);

			  trainsize = nrow(trainingSample);
			  trainsizeCase = nrow(mysampleCases);
			  trainsizeCont = nrow(mysampleControl);
			  
			  sen = sum( 1*((trainingSample [,Outcome] > 0)*( trainmodel$linear.predictor >= 0 )) , na.rm = TRUE);
			  spe = sum( 1*((trainingSample [,Outcome] == 0)*( trainmodel$linear.predictor < 0 )) , na.rm = TRUE);
			  acc = sen+spe;
			  
			  acc = acc/trainsize;
			  sen = sen/trainsizeCase;
			  spe = spe/trainsizeCont;

			  taccuracy <- append(taccuracy,acc);
			  tsensitivity <- append(tsensitivity,sen);
			  tspecificity <- append(tspecificity,spe);

			  cf <- as.vector(trainmodel$coefficients);
			  modelReclas <- getVarReclassification(trainmodel,trainingSample,Outcome=Outcome,type,independentFrame=indpendentSample);

			  
			  wtsIdi <- modelReclas$z.IDIs;
			  for (n in 1:length(modelReclas$z.IDIs))
			  {
					wtsIdi[n] = 0;
					if (!is.na(modelReclas$tz.IDIs[n]) && !is.na(modelReclas$IDIs[n]))
					{
						if (modelReclas$tz.IDIs[n] > 0) 
						{
							wtsIdi[n] <- (modelReclas$IDIs[n]-abs(modelReclas$IDIs[n]-modelReclas$tIDIs[n]))/modelReclas$tIDIs[n];
							if (wtsIdi[n] < -0.1) 
							{
								wtsIdi[n] = -0.1;
							}
						}
					}
			  }
			  if (type != "COX") 
			  {
					wtsIdi <- append(1,wtsIdi)
			  }
			  wtdcf <- wtsIdi*cf;
			  sumwtdcf <- sumwtdcf+wtdcf;
			  sumwts <- sumwts+wtsIdi;
			  valinserted = valinserted +1;
			  
			  if (valinserted == 1) 
			  {
					bcoef <- t(matrix(cf));
					zNRI <- t(matrix(modelReclas$tz.NRIs));
					zIDI <- t(matrix(modelReclas$tz.IDIs));
					test.zNRI <- t(matrix(modelReclas$z.NRIs));
					test.zIDI <- t(matrix(modelReclas$z.IDIs));
					NRI <- t(matrix(modelReclas$NRIs));
					IDI <- t(matrix(modelReclas$IDIs));
			  }
			  else 
			  {
					bcoef <- rbind(bcoef,t(matrix(cf)));	  
					zNRI <- rbind(zNRI,t(matrix(modelReclas$tz.NRIs)));
					zIDI <- rbind(zIDI,t(matrix(modelReclas$tz.IDIs)));
					test.zNRI <- rbind(zNRI,t(matrix(modelReclas$z.NRIs)));
					test.zIDI <- rbind(zIDI,t(matrix(modelReclas$z.IDIs)));
					NRI <- rbind(NRI,t(matrix(modelReclas$NRIs)));
					IDI <- rbind(IDI,t(matrix(modelReclas$IDIs)));

					if (plots && ((doOver %% 20) == 0))
					{				
						par(mfrow=c(1,3))
						hist(taccuracy,main="Accuracy",breaks=20,xlim=c(0.5, 1.0),col="red");
						d <- density(accuracy);  
						lines(d);
						abline(v=gacc/smaptot,col = "green");
						hist(tsensitivity,main="Sensitivity",breaks=20,xlim=c(0.5, 1.0),col="red");
						d <- density(sensitivity)  
						lines(d);
						abline(v=gsen/smapsen,col = "green");
						hist(tspecificity,main="Specificity",breaks=20,xlim=c(0.5, 1.0),col="red");
						d <- density(specificity)  
						lines(d);
						abline(v=gspe/smapspe,col = "green");
					}
			}
		}
		else
		{
			cat(doOver,": Warning Bootstrap fitting error \n");
#			print(trainmodel$coefficients);
		}
		 
		    
	}

	bootmodel <- basemodel;
	bootmodelPond <- basemodel;

	bootmodel$coefficients <- colMedians(bcoef,na.rm = TRUE);

	if (length(sumwts)>0)
	{
		for (n in 1:length(sumwts))
		{
			if (sumwts[n] < 0)  
			{
				sumwts[n] = 1;
				sumwtdcf[n] = 0;  
			}
		}
	}
	
	bootmodelPond$coefficients <- sumwtdcf * (1.0 / sumwts);

	
	names(bootmodel$coefficients) <- names(basemodel$coefficients);
	names(bootmodelPond$coefficients) <- names(basemodel$coefficients);

	pr <- predictForFresa(bootmodel,newdata=dataframe,type = 'linear');
   
	bootmodel$linear.predictors <- pr;

	p <- predictForFresa(bootmodel,newdata=dataframe,type = 'prob');
   
	bootmodel$fitted.values <- p;

	p <- predictForFresa(bootmodelPond,newdata=dataframe,type = 'linear');
   
	bootmodelPond$linear.predictors <- p;

	p <- predictForFresa(bootmodelPond,newdata=dataframe,type = 'prob');
   
	bootmodelPond$fitted.values <- p;
	 
	
	ponacc = 0.0;
	ponsen = 0.0;
	ponspe = 0.0;
	
	sen <- sum( 1*((dataframe[,Outcome] > 0)*( bootmodel$linear.predictors >= 0 )) , na.rm = TRUE)
	spe <- sum( 1*((dataframe[,Outcome] == 0)*( bootmodel$linear.predictors < 0 )) , na.rm = TRUE)
	acc <- sen+spe;

	ponsen <- sum( 1*((dataframe[,Outcome] > 0)*( bootmodelPond$linear.predictors >= 0 )) , na.rm = TRUE)
	ponspe <- sum( 1*((dataframe[,Outcome] == 0)*( bootmodelPond$linear.predictors < 0 )) , na.rm = TRUE)
	ponacc <- ponsen+ponspe;

	acc = acc/framesize;
	sen = sen/sizecases;
	spe = spe/sizecontrol;

	ponacc = ponacc/framesize;
	ponsen = ponsen/sizecases;
	ponspe = ponspe/sizecontrol;

	BlindAccuracy = gacc/smaptot;
	BlindSensitivity = gsen/smapsen;
	BlindSpecificity = gspe/smapspe;
	
	par(mfrow=c(2,2))
	roc(dataframe[,Outcome], bootmodelPond$linear.predictors,col="purple",auc=TRUE,plot=TRUE,smooth=FALSE)
	par(new=TRUE)
	roc( dataframe[,Outcome], basemodel$linear.predictors, col="blue",plot=TRUE,smooth=FALSE);
	par(new=TRUE)
	blidRoc <- roc(testoutcome,testprediction,col="red",auc=TRUE,print.auc=TRUE,plot=TRUE,smooth=FALSE)
	par(new=TRUE)
	bootRoc <- roc( dataframe[,Outcome], bootmodel$linear.predictors,plot=TRUE,ci=plots,auc=TRUE,of='se',specificities=c(0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05),boot.n=200,smooth=FALSE);
	
		plot(ecdf(taccuracy),main="Accuracy",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
		abline(v=BlindAccuracy,col = "red");
		abline(v=acc,col = "blue");
		abline(v=ponacc,col = "green");
		plot(ecdf(tsensitivity),main="Sensitivity",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
		abline(v=BlindSensitivity,col = "red");
		abline(v=sen,col = "blue");
		abline(v=ponsen,col = "green");
		plot(ecdf(tspecificity),main="Specificity",xlim=c(0.5, 1.0),col="black",lwd = 2,verticals = TRUE, do.points = FALSE);
		abline(v=BlindSpecificity,col = "red");
		abline(v=spe,col = "blue");
		abline(v=ponspe,col = "green");

	par(mfrow=c(1,1))
	
	if (class(model.formula) != "formula")
	{
		model.formula = formula(model.formula);
	}
	 
    colnames(IDI) <- attr(terms(model.formula),'term.labels');
    colnames(zIDI) <- attr(terms(model.formula),'term.labels');
    colnames(NRI) <- attr(terms(model.formula),'term.labels');
    colnames(zNRI) <- attr(terms(model.formula),'term.labels');
	colnames(bcoef) <- names(basemodel$coefficients);
	


	result <- structure(llist(
	data=dataframe,
	bin.Outcome=dataframe[,Outcome],
	blind.accuracy=BlindAccuracy,
	blind.sensitivity=BlindSensitivity,
	blind.specificity=BlindSpecificity,
	train.ROCACU=trainRoc,
	blind.ROCAUC= blidRoc,
	boot.ROCAUC=bootRoc, 
	fraction=fraction,
	loops=loops,
	base.Accuracy=baseAcc,
	base.sensitivity=baseSen,
	base.specificity=baseSpe,
	accuracy=accuracy,
	sensitivity=sensitivity,
	specificity=specificity,
	train.accuracy=taccuracy,
	train.sensitivity=tsensitivity,
	train.specificity=tspecificity,
	s.coef=bcoef,
	boot.model=bootmodel,
	mboot.model=bootmodelPond,
	boot.accuracy=acc,
	boot.sensitivity=sen,
	boot.specificity=spe,
	z.NRIs=zNRI,
	z.IDIs=zIDI,
	test.z.NRIs=test.zNRI,
	test.z.IDIs=test.zIDI,
	NRIs=NRI,
	IDIs=IDI,
	blind.test.outcome=testoutcome,
	blind.test.prediction=testprediction,
	labels = FALSE), class = "bootstrapValidation");
	return (result);
}
