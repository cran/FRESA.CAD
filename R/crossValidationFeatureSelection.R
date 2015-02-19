crossValidationFeatureSelection <-
function(size=10,fraction=1.0,pvalue=0.05,loops=100,covariates="1",Outcome,timeOutcome="Time",variableList,data,maxTrainModelSize=10,type=c("LM","LOGIT","COX"),selectionType=c("zIDI","zNRI"),loop.threshold=10,startOffset=0,elimination.bootstrap.steps=25,trainFraction=0.67,trainRepetition=9,elimination.pValue=0.05,CVfolds=10,bootstrap.steps=25,interaction=c(1,1),nk=0,unirank=NULL,print=TRUE,plots=TRUE)
{

if (!requireNamespace("cvTools", quietly = TRUE)) {
   install.packages("cvTools", dependencies = TRUE)
} 

if (!requireNamespace("glmnet", quietly = TRUE)) {
   install.packages("glmnet", dependencies = TRUE)
} 


	enetSamples <- NULL;

	casesample = subset(data,get(Outcome)  == 1);
	controlsample = subset(data,get(Outcome) == 0);

	casesamplesize <- nrow(casesample);

	controlsamplesize <- nrow(controlsample);
	

	K <- as.integer(1.0/(1.0-trainFraction) + 0.5);


	acc = 0.0;
	sen = 0.0;
	spe = 0.0;
	sizecases = 0;
	sizecontrol = 0;
	totsize = 0;
	paracc = 0;
	psen = 0;
	pspe = 0;

	full.acc = 0.0;
	full.sen = 0.0;
	full.spe = 0.0;
	full.paracc = 0;
	full.psen = 0;
	full.pspe = 0;



	formulas <- vector();
	trainCorrelations <- vector();
	blindCorrelations <- vector();
	WholeFoldBlindAccuracy <- vector();
	WholeFoldBlindSpecificity <- vector();
	WholeFoldBlindSensitivity <- vector();
	FoldBlindAccuracy <- vector();
	FoldBlindSpecificity <- vector();
	FoldBlindSensitivity <- vector();
	selection.pValue <- pvalue;
	update.pValue <- c(pvalue,2.0*pvalue);

	CVselection.pValue <- pvalue;
	CVelimination.pValue <- elimination.pValue;
	CVupdate.pValue <- update.pValue


	par(mfrow=c(1,1))
	
	mOrderSel = interaction[1];
	mOrderUpdate = mOrderSel;
	if (length(interaction)>1) mOrderUpdate=interaction[2]
	
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
#	print(shortVarList)
	fullenet <- try(glmnet::cv.glmnet(as.matrix(data[,shortVarList]),as.vector(data[,Outcome]),family="binomial"));
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
	

	
	
	CurModel_Full <- ReclassificationFRESA.Model(size=size,fraction=fraction,pvalue=CVselection.pValue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=data,maxTrainModelSize=maxTrainModelSize,type=type,timeOutcome=timeOutcome,selectionType=selectionType,loop.threshold=loop.threshold,interaction=mOrderSel)
	if (length(CurModel_Full$var.names)==0)
	{
		stop("no initial model found\n");
	}

	if (loops>1)
	{
		UCurModel_Full <- updateModel(Outcome=Outcome,covariates=covariates,pvalue=CVupdate.pValue,VarFrequencyTable=CurModel_Full$ranked.var,variableList=variableList,data=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selectionType,interaction=mOrderUpdate,numberOfModels=1)
	}
	else
	{
		UCurModel_Full <- CurModel_Full;
	}
	if (elimination.bootstrap.steps>1)
	{
		redCurmodel_Full <- bootstrapVarElimination(object=UCurModel_Full$final.model,pvalue=CVelimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=selectionType,loops=elimination.bootstrap.steps,fraction=fraction,print=print,plots=plots);
	}
	else
	{
		
		redCurmodel_Full <- backVarElimination(object=UCurModel_Full$final.model,pvalue=CVelimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=selectionType);
	}
	
	full_formula <- redCurmodel_Full$back.formula;
	fullBootCross <- bootstrapValidation(1.0,bootstrap.steps,full_formula,Outcome,data,type,plots=plots)
	if (is.null(fullBootCross))
	{
		stop("no initial model found\n");
	}
	if (print) summary(fullBootCross,2)

	inserted = 0;
	rocadded = 0;
	split.blindSen <- NULL;
	blindreboot <- NULL;
	kmmSamples <- NULL;
	full.kmmSamples <- NULL;
	totSamples <- NULL;
	full.totSamples <- NULL;

	fullsammples <- min(casesamplesize,controlsamplesize);
	if ( K > fullsammples) K=fullsammples
	cat("Number of folds: ",K,"\n");
	specificities <- c(0.975,0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.05);


	for (i in 1:trainRepetition)
	{
		j <- 1 + ((i-1) %% K)
		if ( j == 1)
		{
#			casefolds <- cvTools::cvFolds(casesamplesize, K, type = "random");
#			controlfolds <- cvTools::cvFolds(controlsamplesize, K, type = "random");
			casefolds <- cvTools::cvFolds(casesamplesize, K,1,  "random");
			controlfolds <- cvTools::cvFolds(controlsamplesize, K,1,  "random");
		}

		CaseTrainSet <- casesample[casefolds$subsets[casefolds$which != j,],];
		CaseBlindSet <- casesample[casefolds$subsets[casefolds$which == j,],];
		ControlTrainSet <- controlsample[controlfolds$subsets[controlfolds$which != j,],];
		ControlBlindSet <- controlsample[controlfolds$subsets[controlfolds$which == j,],];

		TrainSet <- rbind(CaseTrainSet,ControlTrainSet);
		BlindSet <- rbind(CaseBlindSet,ControlBlindSet);
		minTrainSamples <- min(nrow(CaseTrainSet),nrow(ControlTrainSet));

		if (nk==0)
		{
			nk = 2*as.integer(sqrt(minTrainSamples/2)) + 1;
		}
			

		KnnTrainSet <- rbind(CaseTrainSet[sample(1:nrow(CaseTrainSet),minTrainSamples,replace=FALSE),],ControlTrainSet[sample(1:nrow(ControlTrainSet),minTrainSamples,replace=FALSE),])
		
		par(mfrow=c(1,1))

		redBootCross <- bootstrapValidation(1.0,bootstrap.steps,full_formula,Outcome,TrainSet,type,plots=plots)
		par(mfrow=c(1,1))

#		print(summary(redBootCross$boot.model))
		full.p <- predictForFresa(redBootCross$boot.model,BlindSet, 'prob');
		Fullknnclass <- getKNNpredictionFromFormula(redCurmodel_Full$back.formula,KnnTrainSet,BlindSet,Outcome,nk)


		if (!is.null(unirank))
		{
#			cat("Ranking Again \n");
			variableList <- update.uniRankVar(unirank,data=TrainSet,fullAnalysis=FALSE)$orderframe;
			shortVarList <- as.vector(variableList[1:size,1]);
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
		}
		
		
		
		if (!is.null(fullenet))
		{
#			cat("In elastic Net\n")
			foldenet <- try(glmnet::cv.glmnet(as.matrix(TrainSet[,shortVarList]),as.vector(TrainSet[,Outcome]),family="binomial"));
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
#			print(enetVariables)
		}




		cat ("Loop :",i,"Input Cases =",sum(data[,Outcome] > 0 ),"Input Control =",sum(data[,Outcome] == 0),"\n")
		cat ("Loop :",i,"Train Cases =",sum(TrainSet[,Outcome] > 0 ),"Train Control =",sum(TrainSet[,Outcome] == 0),"\n")
		cat ("Loop :",i,"Blind Cases =",sum(BlindSet[,Outcome] > 0 ),"Blind Control =",sum(BlindSet[,Outcome] == 0),"\n")
		cat ("K   :",nk,"KNN T Cases =",sum(KnnTrainSet[,Outcome] > 0 ),"KNN T Control =",sum(KnnTrainSet[,Outcome] == 0),"\n")

		CurModel_S <- ReclassificationFRESA.Model(size=size,fraction=fraction,pvalue=selection.pValue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=variableList,data=TrainSet,maxTrainModelSize=maxTrainModelSize,type=type,timeOutcome=timeOutcome,selectionType=selectionType,loop.threshold=loop.threshold,interaction=mOrderSel)
		if (length(CurModel_S$var.names)>0)
		{
			if (loops>1)
			{
				UCurModel_S <- updateModel(Outcome=Outcome,covariates=covariates,pvalue=update.pValue,VarFrequencyTable=CurModel_S$ranked.var,variableList=variableList,data=TrainSet,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selectionType,interaction=mOrderUpdate,numberOfModels=1)
			}
			else
			{
				UCurModel_S <- CurModel_S;
			}
			if (elimination.bootstrap.steps > 1)
			{
				redCurmodel_S <- bootstrapVarElimination(object=UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=TrainSet,startOffset=startOffset,type=type,selectionType=selectionType,loops=elimination.bootstrap.steps,fraction=fraction,print=print,plots=plots);	
				redBootCross_S <- redCurmodel_S$bootCV;
			}
			else
			{
				redCurmodel_S <- backVarElimination(object=UCurModel_S$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=TrainSet,startOffset=startOffset,type=type,selectionType=selectionType);
				redBootCross_S <- bootstrapValidation(1.0,bootstrap.steps,redCurmodel_S$back.formula,Outcome,TrainSet,type,plots=plots)
				par(mfrow=c(1,1))
			}
			 

			cat ("\n The last CV bootstrapped model")
			if (print) s <- summary(redBootCross_S,2)

			if (redCurmodel_S$lastRemoved >= 0)
			{
#				print(summary(redCurmodel_S$back.model))
				p <- predictForFresa(redBootCross_S$boot.model,BlindSet, 'prob');
				knnclass <- getKNNpredictionFromFormula(redCurmodel_S$back.formula,KnnTrainSet,BlindSet,Outcome,nk)
				
								
				inserted = inserted + 1

				tcor <- cor.test(p, full.p, method = "spearman",na.action=na.omit,exact=FALSE)$estimate
				trainCorrelations <- append(trainCorrelations,tcor);
				framesize <- nrow(BlindSet);
				bcor <- 0;
				if (framesize>5)
				{
					bcor <- cor.test(p, full.p, method = "spearman",na.action=na.omit,exact=FALSE)$estimate;
					blindCorrelations <- append(blindCorrelations,bcor);
				}

				if ((sum(BlindSet[,Outcome]>0) > 3) && (sum(BlindSet[,Outcome]==0) > 3))
				{
					splitRoc <- pROC::roc(BlindSet[,Outcome], p,plot=FALSE,ci=TRUE,auc=TRUE,of='se',specificities=specificities,boot.n=100,smooth=FALSE)$ci[,2]
					if (rocadded == 0)
					{
						split.blindSen <- splitRoc;
					}
					else
					{
						split.blindSen <- rbind(split.blindSen,splitRoc);
					}
					rocadded = rocadded + 1;
				}
				
				totsize <- totsize + framesize;
				scase <- sum(BlindSet[,Outcome] == 1);
				scotr <- sum(BlindSet[,Outcome] == 0);
				sizecases <- sizecases + scase;
				sizecontrol <- sizecontrol + scotr;
				psen <- sum( 1*((BlindSet[,Outcome] > 0)*( p >= 0.5 )) , na.rm = TRUE)
				pspe <- sum( 1*((BlindSet[,Outcome] == 0)*( p < 0.5 )) , na.rm = TRUE)
				acc <- acc + psen + pspe;
				sen <- sen + psen;
				spe <- spe + pspe;
				psen <- sum( 1*((BlindSet[,Outcome] > 0)*( full.p >= 0.5 )) , na.rm = TRUE)
				pspe <- sum( 1*((BlindSet[,Outcome] == 0)*( full.p < 0.5 )) , na.rm = TRUE)
				full.acc <- full.acc + psen + pspe;
				full.sen <- full.sen + psen;
				full.spe <- full.spe + pspe;
				paracc = acc/totsize;
				psen = 0;
				pspe = 0;
				if (sizecases>0) 
				{
					psen = sen/sizecases;
				}
				if (sizecontrol>0) 
				{
					pspe = spe/sizecontrol;
				}

				full.paracc = full.acc/totsize;
				full.psen = 0;
				full.pspe = 0;
				if (sizecases>0) 
				{
					full.psen = full.sen/sizecases;
				}
				if (sizecontrol>0) 
				{
					full.pspe = full.spe/sizecontrol;
				}

				WholeFoldBlindAccuracy <- append(WholeFoldBlindAccuracy,redBootCross$blind.accuracy);
				WholeFoldBlindSpecificity <- append(WholeFoldBlindSpecificity,redBootCross$blind.specificity);
				WholeFoldBlindSensitivity <- append(WholeFoldBlindSensitivity,redBootCross$blind.sensitivity);

				FoldBlindAccuracy <- append(FoldBlindAccuracy,redBootCross_S$blind.accuracy);
				FoldBlindSpecificity <- append(FoldBlindSpecificity,redBootCross_S$blind.specificty);
				FoldBlindSensitivity <- append(FoldBlindSensitivity,redBootCross_S$blind.sensitivity);

				if (inserted == 1)
				{
					totSamples <- cbind(BlindSet[,Outcome],p,i);
					full.totSamples <- cbind(BlindSet[,Outcome],full.p,i);
					kmmSamples <- cbind(BlindSet[,Outcome],abs(knnclass$prob$prob-1*(knnclass$prediction=="0")),i);
					full.kmmSamples <- cbind(BlindSet[,Outcome],abs(Fullknnclass$prob$prob-1*(Fullknnclass$prediction=="0")),i);
				}
				else
				{
					px <- cbind(BlindSet[,Outcome],p,i);
					totSamples <- rbind(totSamples,px);
					px <- cbind(BlindSet[,Outcome],full.p,i);
					full.totSamples <- rbind(full.totSamples,px);
					px <- cbind(BlindSet[,Outcome],abs(knnclass$prob$prob-1*(knnclass$prediction=="0")),i);
					kmmSamples <- rbind(kmmSamples,px);
					px <- cbind(BlindSet[,Outcome],abs(Fullknnclass$prob$prob-1*(Fullknnclass$prediction=="0")),i);
					full.kmmSamples <- rbind(full.kmmSamples,px);
				}
				
				formulas <- append(formulas,as.character(redCurmodel_S$back.formula));
				knnACC <- sum(kmmSamples[,1] == (kmmSamples[,2]>0.5))/totsize;
				knnSEN <- sum((kmmSamples[,1]>0.5) & (kmmSamples[,2]>0.5))/sizecases;
				knnSPE <- sum((kmmSamples[,1]<0.5) & (kmmSamples[,2]<0.5))/sizecontrol;

				full.knnACC <- sum(full.kmmSamples[,1] == (full.kmmSamples[,2]>0.5))/totsize;
				full.knnSEN <- sum((full.kmmSamples[,1]>0.5) & (full.kmmSamples[,2]>0.5))/sizecases;
				full.knnSPE <- sum((full.kmmSamples[,1]<0.5) & (full.kmmSamples[,2]<0.5))/sizecontrol;


				cat ("Loop :",i,"Blind Cases =",scase,"Blind Control =",scotr,"Total =",totsize, "Size Cases =",sizecases,"Size Control =",sizecontrol,"\n")
				cat ("Accumulated Models CV Accuracy        =",paracc,"Sensitivity =",psen,"Specificity =",pspe,"\n")
				cat ("Initial Model Accumulated CV Accuracy =",full.paracc,"Sensitivity =",full.psen,"Specificity =",full.pspe,"\n");
				cat ("Initial Model Bootstrapped Accuracy   =",redBootCross$blind.accuracy,"Sensitivity =",redBootCross$blind.sensitivity,"Specificity =",redBootCross$blind.specificity,"\n")
				cat ("Current Model Bootstrapped Accuracy   =",redBootCross_S$blind.accuracy,"Sensitivity =",redBootCross_S$blind.sensitivity,"Specificity =",redBootCross_S$blind.specificity,"\n")
				cat ("Current KNN Accuracy   =",knnACC,"Sensitivity =",knnSEN,"Specificity =",knnSPE,"\n")
				cat ("Initial KNN Accuracy   =",full.knnACC,"Sensitivity =",full.knnSEN,"Specificity =",full.knnSPE,"\n")
				cat ("Train Correlation: ",tcor," Blind Correlation :",bcor,"\n KNN to Model Confusion Matrix: \n")
				print(table(kmmSamples[,2]>0.5,totSamples[,2]>0.5))
			}
			else
			{
				cat ("Loop :",i,"No Model.\n")
			}
		}
		else
		{
			cat ("Loop :",i,"No Model.\n")
		}
		
	}
	if (length(formulas)==0)
	{
		stop("No Significant Models Found\n");
	}
	colnames(totSamples) <- c("Outcome","Blind_Prediction","Model");
	totSamples <- as.data.frame(totSamples);
	colnames(full.totSamples) <- c("Outcome","Blind_Prediction","Model");
	full.totSamples <- as.data.frame(full.totSamples);

	colnames(kmmSamples) <- c("Outcome","Blind_Prediction","Model");
	kmmSamples <- as.data.frame(kmmSamples);
	
	colnames(full.kmmSamples) <- c("Outcome","Blind_Prediction","Model");
	full.kmmSamples <- as.data.frame(full.kmmSamples);

	if (!is.null(enetSamples))
	{
		colnames(enetSamples) <- c("Outcome","Blind_Prediction","Model");
		enetSamples <- as.data.frame(enetSamples);
	}

	
	
	plotModels.ROC(totSamples);
	par(mfrow=c(1,1))
	incBsen=0
	aucBlindTest <- pROC::roc(totSamples[,1],totSamples[,2],col="red",auc=TRUE,plot=TRUE,smooth=FALSE,lty=3)$auc
	par(new=TRUE)
	aucCVBlind <- pROC::roc(full.totSamples[,1],full.totSamples[,2],col="blue",auc=TRUE,plot=TRUE,ci=TRUE,smooth=FALSE)$auc
	par(new=TRUE)
	aucTrain <- pROC::roc( fullBootCross$outcome, fullBootCross$boot.model$linear.predictors,col="green",plot=TRUE,auc=TRUE,smooth=FALSE)$auc;        
	par(new=TRUE)
	aucBoot <- pROC::roc( fullBootCross$testOutcome, fullBootCross$testPrediction,col="black",auc=TRUE,plot=TRUE,smooth=FALSE)$auc;
	ley.names <- c(paste("Bootstrapped: Train Model ROC (",sprintf("%.3f",aucTrain),")"),paste("Bootstrapped: Blind ROC (",sprintf("%.3f",aucBoot),")"),
	paste("CV: Blind ROC (",sprintf("%.3f",aucCVBlind),")"),paste("CV: Blind Fold Models Coherence (",sprintf("%.3f",aucBlindTest),")"))
	ley.colors <- c("green","black","blue","red")
	ley.lty <- c(1,1,1,3)
	if (rocadded>0)
	{
		boxplot(split.blindSen,add=TRUE, axes = FALSE,boxwex=0.04,at=specificities);
		sumSen <- colMeans(split.blindSen,na.rm = TRUE);
		sennames <- names(sumSen);
		sumSen <- append(0,sumSen);
		sumSen <- append(sumSen,1);
		sennames <- append("1",sennames);
		sennames <- append(sennames,"0");
		names(sumSen) <- sennames;
		spevalues <- as.numeric(names(sumSen));
		lines(spevalues,sumSen,col="red",lwd=2.0);
		auc = 0;
		for (i in 2:length(spevalues))
		{
			auc = auc + (spevalues[i-1]-spevalues[i])*(sumSen[i-1]+(sumSen[i]-sumSen[i-1])/2)
		}
		ley.names <- append(ley.names,paste("CV Blind: Mean ROC of Models (",sprintf("%.3f",auc),")"));
		ley.colors <- append(ley.colors,"red");
		ley.lty  <- append(ley.lty,1);
	}
	else
	{
		sumSen = NA;
	}
	
	legend(0.6,0.30, legend=ley.names,col = ley.colors, lty = ley.lty,bty="n")
	
	result <- list(formula.list=formulas,
	Models.testPrediction=totSamples,
	FullModel.testPrediction=full.totSamples,
	TestRetrained.blindPredictions=blindreboot,
	LastTrainedModel.bootstrapped=redCurmodel_S$bootCV,
	Test.accuracy=paracc,
	Test.sensitivity=psen,
	Test.specificity=pspe,
	Train.correlationsToFull=trainCorrelations,
	Blind.correlationsToFull=blindCorrelations,
	FullModelAtFoldAccuracies=WholeFoldBlindAccuracy,
	FullModelAtFoldSpecificties=WholeFoldBlindSpecificity,
	FullModelAtFoldSensitivities=WholeFoldBlindSensitivity,
	AtCVFoldModelBlindAccuracies=FoldBlindAccuracy,
	AtCVFoldModelBlindSpecificities=FoldBlindSpecificity,
	AtCVFoldModelBlindSensitivities=FoldBlindSensitivity,
	Models.CVblindMeanSensitivites=sumSen,
	varIDISelection = CurModel_Full,
	updateIDISelection = UCurModel_Full,
	backIDIElimination = redCurmodel_Full,
	FullModel.bootstrapped=fullBootCross,
	Models.testSensitivities = split.blindSen,
	FullKNN.testPrediction=full.kmmSamples,
	KNN.testPrediction=kmmSamples,
	fullenet=fullenet,
	enet.testPredictions=enetSamples,
	enetVariables=enetVariables);
	return (result)
}
