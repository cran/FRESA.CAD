crossValidationNeRIFeatureSelection <-
function(size=10,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,timeOutcome="Time",variableList,data,maxTrainModelSize=10,type=c("LM","LOGIT","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loop.threshold=10,startOffset=0,elimination.bootstrap.steps=25,trainFraction=0.67,trainRepetition=9,elimination.pValue=0.05,setIntersect=1,interaction=c(1,1),update.pvalue=c(0.05,0.05),unirank=NULL,print=TRUE,plots=TRUE)
{

if (!requireNamespace("cvTools", quietly = TRUE)) {
   install.packages("cvTools", dependencies = TRUE)
} 

if (!requireNamespace("glmnet", quietly = TRUE)) {
   install.packages("glmnet", dependencies = TRUE)
} 


	enetSamples <- NULL;
	K <- as.integer(1.0/(1.0-trainFraction) + 0.5);
	mOrderSel = interaction[1];
	mOrderUpdate = mOrderSel;
	if (length(interaction)>1) 
	{
		mOrderUpdate=interaction[2];
	}

	cat(type,"\n")

	fullsammples <- nrow(data);
	if ( K > fullsammples) K=fullsammples


	shortVarList <- as.vector(variableList[1:size,1]);
	varlist <- vector();
	for (i in 1:length(shortVarList))
	{
		varlist <- append(varlist,str_replace_all(unlist(strsplit(
						str_replace_all(
							str_replace_all(
								str_replace_all(
									str_replace_all(shortVarList[i],"I\\("," ")
								,"\\("," ")
							,">","\\*")
						,"<","\\*")
				,"\\*"))[1]," ",""))
	}
	shortVarList <- as.vector(rownames(table(varlist)))
	if (type=="LM")
	{
		fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="gaussian"));
	}
	else
	{
		fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="binomial"));
	}
	if (inherits(fullenet, "try-error"))
	{
		cat("enet Error")
		fullenet <- NULL;
	}
	else
	{
		cenet <- as.matrix(coef(fullenet))
		print(enetVariables <- list(names(cenet[as.vector(cenet[,1]>0),])))
	}


	Full_CurModel_S <- NeRIBasedFRESA.Model(size=size,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,loop.threshold=loop.threshold,interaction=mOrderSel)
	Full_UCurModel_S <- NULL;
	Full_redCurmodel_S <- NULL;
	if (length(Full_CurModel_S$var.names)==0)
	{
		stop("no model found\n");
	}
	if (loops>1)
	{
		Full_UCurModel_S <- updateNeRIModel(Outcome=Outcome,covariates=covariates,pvalue=update.pvalue,VarFrequencyTable=Full_CurModel_S$ranked.var,variableList=variableList,data=data,type=type,testType=testType,timeOutcome=timeOutcome,interaction=mOrderUpdate)
	}
	else
	{
		Full_UCurModel_S <- Full_CurModel_S;
	}
	if (elimination.bootstrap.steps == 1 )
	{
		Full_redCurmodel_S <- backVarNeRIElimination(object=Full_UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,setIntersect=setIntersect);
		Full_redCurmodel_S$bootCV  <- bootstrapValidationNeRI(1,100,Full_redCurmodel_S$back.formula,Outcome,data,type,plots=plots)
			
	}
	else
	{		
#		cat("In Reduction \n")
		Full_redCurmodel_S <- bootstrapVarNeRIElimination(object=Full_UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,
													data=data,startOffset=startOffset,type=type,
													testType=testType,loops=elimination.bootstrap.steps,setIntersect=setIntersect,print=print,plots=plots);
	}

#	cat(format(Full_redCurmodel_S$back.formula)," <-Back formula\n");
#	print(summary(Full_redCurmodel_S$back.model));
	
	
	formulas <- vector();

	vtrainRMS <- vector();
	vblindRMS <- vector();

	vtrainSpearman <- vector();
	vtrainPearson <- vector();

	FullvtrainRMS <- vector();
	FullvblindRMS <- vector();
	FullvtrainSpearman <- vector();
	FullvtrainPearson <- vector();

	blindFoldPearson <-  vector();
	blindFoldSpearman <- vector();
	blindFoldCstat <- vector();
	blindFoldMS <- vector();
	CVBlindPearson  <- vector();
	CVBlindSpearman  <- vector();
	CVBlindRMS <- vector();

	inserted = 0;
	lastmodel = 0;

	for (i in 1:trainRepetition)
	{
		j <- 1 + ((i-1) %% K)
		if ( j == 1)
		{
#			lowFolds <- cvTools::cvFolds(nrow(lowsample), K, type = "random");
#			highFolds <- cvTools::cvFolds(nrow(hihgsample), K, type = "random");
#			lowFolds <- cvTools::cvFolds(nrow(lowsample), K,1,  "random");
#			highFolds <- cvTools::cvFolds(nrow(hihgsample), K,1, "random");
			sampleFolds <- cvTools::cvFolds(nrow(data), K,1, "random");
		}
		
#		lowTrainSet <- lowsample[lowFolds$subsets[lowFolds$which != j,],];
#		highTainSet <- hihgsample[highFolds$subsets[highFolds$which != j,],];
#		lowTestSet <- lowsample[lowFolds$subsets[lowFolds$which == j,],];
#		highTestSet <- hihgsample[highFolds$subsets[highFolds$which == j,],];
		
#		TrainSet <- rbind(lowTrainSet,highTainSet);
#		BlindSet <- rbind(lowTestSet,highTestSet);
		

		TrainSet <- data[sampleFolds$subsets[sampleFolds$which != j,],];
		BlindSet <- data[sampleFolds$subsets[sampleFolds$which == j,],];

		if (!is.null(unirank))
		{
			variableList <- update.uniRankVar(unirank,data=TrainSet,fullAnalysis=FALSE)$orderframe;
			shortVarList <- as.vector(variableList[1:size,1]);
#			print(shortVarList)
			varlist <- vector();
			for (nn in 1:length(shortVarList))
			{
				varlist <- append(varlist,str_replace_all(unlist(strsplit(
								str_replace_all(
									str_replace_all(
										str_replace_all(
											str_replace_all(shortVarList[nn],"I\\("," ")
										,"\\("," ")
									,">","\\*")
								,"<","\\*")
						,"\\*"))[1]," ",""))
			}
			shortVarList <- as.vector(rownames(table(varlist)))
#			print(shortVarList)
		}




		cat("Samples Train :",nrow(TrainSet),"Samples Test :",nrow(BlindSet),"\n");
		cat ("Loop :",i,"\n")

		if (!is.null(fullenet))
		{
			if (type=="LM")
			{
				foldenet <- try(glmnet::cv.glmnet(as.matrix(TrainSet[,shortVarList]),as.vector(TrainSet[,Outcome]),family="gaussian"));
			}
			else
			{
				foldenet <- try(glmnet::cv.glmnet(as.matrix(TrainSet[,shortVarList]),as.vector(TrainSet[,Outcome]),family="binomial"));
			}
			cenet <- as.matrix(coef(foldenet))
			enetVariables[[i+1]] <- names(cenet[as.vector(cenet[,1]>0),])
			if (i == 1)
			{
				enetSamples <- cbind(BlindSet[,Outcome],predict(foldenet,as.matrix(BlindSet[,shortVarList])),i);
			}
			else
			{
				enetSamples <- rbind(enetSamples,cbind(BlindSet[,Outcome],predict(foldenet,as.matrix(BlindSet[,shortVarList])),i));
			}
##			print(enetVariables)
		}

		
		
		par(mfrow=c(1,1))

		cat(type,"\n")

#		cat(Full_redCurmodel_S$back.formula," <-Back formula\n");

		CurModel_S <- NeRIBasedFRESA.Model(size=size,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=TrainSet,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,loop.threshold=loop.threshold,interaction=mOrderSel)
		if (length(CurModel_S$var.names)>0)
		{
			inserted = inserted +1;
			if (loops>1)
			{
				UCurModel_S <- updateNeRIModel(Outcome=Outcome,covariates=covariates,pvalue=update.pvalue,VarFrequencyTable=CurModel_S$ranked.var,variableList=variableList,data=TrainSet,type=type,testType=testType,timeOutcome=timeOutcome,interaction=mOrderUpdate)
			}
			else
			{
				UCurModel_S <- CurModel_S;
			}
			if (elimination.bootstrap.steps == 1 )
			{
				redCurmodel_S <- backVarNeRIElimination(object=UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=TrainSet,startOffset=startOffset,type=type,testType=testType,setIntersect=setIntersect);
				redCurmodel_S$bootCV  <- bootstrapValidationNeRI(1,100,redCurmodel_S$back.formula,Outcome,TrainSet,type,plots=plots)
			}
			else
			{
				redCurmodel_S <- bootstrapVarNeRIElimination(object=UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=TrainSet,startOffset=startOffset,type=type,testType=testType,loops=elimination.bootstrap.steps,setIntersect=setIntersect,print=print,plots=plots);
			}

#			cat(Full_redCurmodel_S$back.formula," -> After Model\n");
			Full_model <- modelFitting(Full_redCurmodel_S$back.formula,TrainSet,type)	
#			print(summary(Full_model))
#			cat("After refitting\n");
			
			predictTest <- predictForFresa(redCurmodel_S$bootCV$boot.model,BlindSet,"linear");
			Full_predictTest <- predictForFresa(Full_model,BlindSet,"linear");

#						cat("After Prediction\n");

			trainResiduals <- residualForNeRIs(redCurmodel_S$bootCV$boot.model,TrainSet,Outcome);
			blindResiduals <- residualForNeRIs(redCurmodel_S$bootCV$boot.model,BlindSet,Outcome);
			FulltrainResiduals <- residualForNeRIs(Full_model,TrainSet,Outcome);
			FullblindResiduals <- residualForNeRIs(Full_model,BlindSet,Outcome);

#									cat("After Residuals\n");

			if (inserted == 1)
			{
				totSamples <- cbind(BlindSet[,Outcome],predictTest,i,blindResiduals);
				full.totSamples <- cbind(BlindSet[,Outcome],Full_predictTest,i,FullblindResiduals);
			}
			else
			{
				px <- cbind(BlindSet[,Outcome],predictTest,i,blindResiduals);
				totSamples <- rbind(totSamples,px);
				px <- cbind(BlindSet[,Outcome],Full_predictTest,i,FullblindResiduals);
				full.totSamples <- rbind(full.totSamples,px);
			}
			formulas <- append(formulas,as.character(redCurmodel_S$back.formula));

			trainPredictTest <- predictForFresa(redCurmodel_S$bootCV$boot.model,TrainSet,"linear");
			trainFull_predictTest <- predictForFresa(Full_model,TrainSet,"linear");
				
			
			trainRMS <- sqrt(sum(trainResiduals^2)/nrow(TrainSet));
			trainPearson <- cor.test(TrainSet[,Outcome], trainPredictTest, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
			trainSpearman <- cor.test(TrainSet[,Outcome], trainPredictTest, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

			FulltrainRMS <- sqrt(sum(FulltrainResiduals^2)/nrow(TrainSet));
			FulltrainPearson <- cor.test(TrainSet[,Outcome], trainFull_predictTest, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
			FulltrainSpearman <- cor.test(TrainSet[,Outcome], trainFull_predictTest, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

			blindRMS <- sqrt(sum(blindResiduals^2)/nrow(BlindSet));
			FullblindRMS <- sqrt(sum(FullblindResiduals^2)/nrow(BlindSet));

			if (nrow(BlindSet)>5)
			{
				foldPearson <- cor.test(BlindSet[,Outcome], predictTest, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
				foldSpearman <- cor.test(BlindSet[,Outcome], predictTest, method = "spearman",na.action=na.omit,exact=FALSE)$estimate				
				cstat <- rcorr.cens(predictTest,BlindSet[,Outcome], outx=FALSE)[1];
				foldRMS <- sum(blindResiduals^2)/(nrow(BlindSet)-1);
				cat("Fold RMS: ",sqrt(foldRMS),"Fold Test Pearson: ", foldPearson, "Fold Test Spearman: ",foldSpearman,"Fold Cstat:,",cstat,"\n");

				blindFoldMS <- append(blindFoldMS,foldRMS);
				blindFoldPearson <- append(blindFoldPearson,foldPearson);
				blindFoldSpearman <- append(blindFoldSpearman,foldSpearman);
				blindFoldCstat <- append(blindFoldCstat,cstat);
				cat("Accu RMS: ",sqrt(mean(blindFoldMS)),"Accu Test Pearson: ", mean(blindFoldPearson), "Accu Test Spearman: ",mean(blindFoldSpearman),"Accu Cstat:,",mean(blindFoldCstat),"\n");
			}
#			cat("After RMS\n");

			if (nrow(totSamples)>5)
			{
				blindPearson <- cor.test(totSamples[,1], totSamples[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate
				blindSpearman <- cor.test(totSamples[,1], totSamples[,2], method = "spearman",na.action=na.omit,exact=FALSE)

				cstat <- rcorr.cens(totSamples[,2],totSamples[,1], outx=FALSE)[1];

				FullblindPearson <- cor.test(full.totSamples[,1], full.totSamples[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate
				FullblindSpearman <- cor.test(full.totSamples[,1], full.totSamples[,2], method = "spearman",na.action=na.omit,exact=FALSE)

				AcumRMS <- sqrt(sum(totSamples[,4]^2)/nrow(totSamples));
				AcumFullRMS <- sqrt(sum(full.totSamples[,4]^2)/nrow(totSamples));


				cat("Samples: ", nrow(totSamples),"Full RMS:",AcumFullRMS," Accumulated Blind RMS: ", AcumRMS," c-index : ",cstat,"\n");
				cat("Full Blind RMS: ", FullblindRMS, " Full Train RMS: ",FulltrainRMS,"\n");
				cat("Blind Pearson: ", blindPearson, " Train Pearson: ",trainPearson,"\n");
				cat("Blind Spearman: ", blindSpearman$estimate,"(", blindSpearman$p.value,") Train Spearman: ",trainSpearman,"\n");
				cat("Full Blind Pearson: ", FullblindPearson , " Full Train Pearson: ",FulltrainPearson,"\n");
				cat("Full Blind Spearman: ", FullblindSpearman$estimate, "(",FullblindSpearman$p.value,") Full Train Spearman: ",FulltrainSpearman,"\n");
			}
			
			vblindRMS <- append(vblindRMS,blindRMS);
			FullvblindRMS <- append(FullvblindRMS,FullblindRMS);
			
			vtrainRMS <- append(vtrainRMS,trainRMS);
			vtrainSpearman <- append(vtrainSpearman,trainSpearman);
			vtrainPearson <- append(vtrainPearson,trainPearson);

			FullvtrainRMS <- append(FullvtrainRMS,FulltrainRMS);
			FullvtrainSpearman <- append(FullvtrainSpearman,FulltrainSpearman);
			FullvtrainPearson <- append(FullvtrainPearson,FulltrainPearson);


		}
		if ( (i %% K) == 0)
		{
			foldtest <- totSamples[totSamples[,3]>lastmodel,];
			CVBlindPearson <- append(CVBlindPearson,cor.test(foldtest[,1], foldtest[,2], method = "pearson",na.action=na.omit,exact=FALSE)$estimate);
			CVBlindSpearman <- append(CVBlindSpearman,cor.test(foldtest[,1], foldtest[,2], method = "spearman",na.action=na.omit,exact=FALSE)$estimate);
			CVBlindRMS <- append(CVBlindRMS,sqrt(sum((foldtest[,1]-foldtest[,2])^2)/nrow(foldtest)));
			lastmodel = i;
		}
	}
			
#	print(enetVariables)

	colnames(totSamples) <- c("Outcome","Blind_Prediction","Model","Residuals");
	totSamples <- as.data.frame(totSamples);

	colnames(full.totSamples) <- c("Outcome","Blind_Prediction","Model","Residuals");
	full.totSamples <- as.data.frame(full.totSamples);

	if (!is.null(enetSamples))
	{
		colnames(enetSamples) <- c("Outcome","Blind_Prediction","Model");
		enetSamples <- as.data.frame(enetSamples);
	}

	blindRMS <- sqrt(sum((totSamples$Residuals)^2)/nrow(totSamples));
	blindPearson <- cor.test(totSamples$Outcome, totSamples$Blind_Prediction, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
	blindSpearman <- cor.test(totSamples$Outcome, totSamples$Blind_Prediction, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

	FullblindRMS <- sqrt(sum((full.totSamples$Residuals)^2)/nrow(totSamples));
	FullblindPearson <- cor.test(full.totSamples$Outcome, full.totSamples$Blind_Prediction, method = "pearson",na.action=na.omit,exact=FALSE)$estimate
	FullblindSpearman <- cor.test(full.totSamples$Outcome, full.totSamples$Blind_Prediction, method = "spearman",na.action=na.omit,exact=FALSE)$estimate

	cstat <- rcorr.cens(totSamples$Blind_Prediction,totSamples$Outcome, outx=FALSE)[1];

	cat("##### CV of Initial Model (Biased) ###### \n");
	cat("Full Blind RMS: ", FullblindRMS,"\n")
	cat("Full Blind Spearman: ", FullblindSpearman,"\n")
	cat("Full Blind Pearson: ", FullblindPearson,"\n")
	if (length(blindFoldMS)>0)
	{
		cat("##### By Fold Analysis ###### \n Samples: ",length(blindFoldMS),"\n");
		cat("Sampled RMS: ", sqrt(mean(blindFoldMS)),"\n")
		cat("Mean Blind Spearman: ", mean(blindFoldSpearman)," (",sd(blindFoldSpearman),")\n")
		cat("Mean Blind Pearson: ", mean(blindFoldPearson)," (",sd(blindFoldPearson),")\n")
		cat("Mean Blind cstat: ",mean(blindFoldCstat)," (",sd(blindFoldCstat),")\n")
	}
	cat("##### Full Set Coherence Analysis ###### \n");
	cat("Blind  RMS: ", blindRMS," c-index : ",cstat,"\n");
	cat("Blind  Spearman: ", blindSpearman,"\n")
	cat("Blind  Pearson: ", blindPearson,"\n")
	if (length(CVBlindPearson)>0)
	{
		cat("##### Models Coherence Analysis ###### \n Samples: ",length(CVBlindPearson),"\n");
		cat("Mean Blind Spearman: ", mean(CVBlindSpearman)," (",sd(CVBlindSpearman),")\n")
		cat("Mean Blind Pearson: ", mean(CVBlindPearson)," (",sd(CVBlindPearson),")\n")
	}

	
	
	result <- list(		formula.list=formulas, 
						Models.testPrediction=totSamples,
						FullModel.testPrediction=full.totSamples,
						backNeRIElimination=Full_redCurmodel_S,
						varNeRISelection=Full_CurModel_S,
						updateNeRISelection=Full_UCurModel_S,
						testRMSE = blindRMS,
						testPearson = blindPearson,
						testSpearman = blindSpearman,
						fullTestRMSE = FullblindRMS,
						fullTestPearson = FullblindPearson,
						fullTestSpearman = FullblindSpearman,
						trainRMSE = vtrainRMS,
						trainPearson = vtrainPearson,
						trainSpearman = vtrainSpearman,
						fullTrainRMS = FullvtrainRMS,
						fullTrainPearson = FullvtrainPearson,
						fullTrainSpearman = FullvtrainSpearman,
						byFoldTestMS = blindFoldMS,
						byFoldTestSpearman = blindFoldSpearman,
						byFoldTestPearson = blindFoldPearson,
						byFoldCstat = blindFoldCstat,
						testRMSEAtFold = vblindRMS,
						FullTestRMSEAtFold = FullvblindRMS,
						fullenet=fullenet,
						enet.testPredictions=enetSamples,
						enetVariables=enetVariables,
						CVBlindPearson=CVBlindPearson,
						CVBlindSpearman=CVBlindSpearman,
						CVBlindRMS=CVBlindRMS
					);
	return (result)
}
