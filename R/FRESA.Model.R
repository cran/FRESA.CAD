FRESA.Model <-
function(formula,data,OptType=c("Binary","Residual"),pvalue=0.05,filter.p.value=0.10,loops=1,maxTrainModelSize=10,loop.threshold=20,elimination.bootstrap.steps=100,bootstrap.steps=100,interaction=c(1,1),print=TRUE,plot=TRUE,folds=10,repeats=1,k=0,categorizationType=c("Raw","Categorical","ZCategorical","RawZCategorical","RawZTail","RawTail"),cateGroups=c(0.1,0.9),raw.dataFrame=NULL,var.description=NULL,testType=c("zIDI","zNRI","Binomial","Wilcox","tStudent","Ftest"))
{
	
	relaxedUpdatePvalue = 2.0;
	cvObject <-  NULL;
	univariate <- NULL;
	if (class(formula)=="character")
	{
		formula <- formula(formula);
	}
	if (class(formula)=="formula")
	{
		OptType <- match.arg(OptType)
		
		varlist <- attr(terms(formula),"variables")
		
		dependent <- as.character(varlist[[2]])
		
		timeOutcome = NA;
		Outcome = NA;
		
		type = "LOGIT";
		if (length(dependent)==3)
		{
			type = "COX"
			timeOutcome = dependent[2];
			Outcome = dependent[3];
			dependentout = paste(dependent[1],"(",dependent[2],",",dependent[3],")");
		}
		else
		{
			Outcome = dependent[1];
			dependentout = Outcome;	
		}
		
		setIntersect <- attr(terms(formula),"intercept")
		if (setIntersect == 0)
		{
			covariates = "0";
		}
		else
		{
			covariates = "1";
		}
		
		termslist <- attr(terms(formula),"term.labels");
		cat (termslist,"\n")
		if (length(termslist)>0)
		{
			for (i in 1:length(termslist))
			{
				covariates <- paste(covariates," + ",termslist[i]);
			}
		}
		cat (covariates,"\n")

		startOffset = length(termslist);
		variables <- vector();
		descrip <- vector();
		pnames <- as.vector(colnames(data));
		for (i in 1:length(pnames))
		{
			detected = 0;
			if (length(termslist)>0)
			{
				for (j in 1:length(termslist))
				{
					if (termslist[j] == pnames[i]) detected = 1;
				}
			}
			if (Outcome == pnames[i]) detected = 1;
			if (!is.na(timeOutcome) )
			{
				if (timeOutcome == pnames[i]) detected = 1;
			}
			if (detected == 0)
			{
				variables <- append(variables,pnames[i]);
				if (!is.null(var.description))
				{
					descrip <- append(descrip,var.description[i]);
				}
			}
		}
		if (!is.null(var.description))
		{
			variables <- cbind(variables,descrip);
		}
		else
		{
			variables <- cbind(variables,variables);
		}
		
		colnames(variables) <- c("Var","Description");
	
		
		trainFraction <- 1.0-1.0/folds;
		trainRepetition <- repeats*folds;
		fraction = 1.0;
		varMax = nrow(variables);
		baseModel <- paste(dependentout,"~",covariates);
		cvObject = NULL;
		reducedModel = NULL;
		if (OptType == "Binary")
		{

			if (length(dependent)==1)
			{
				type = "LOGIT";
			}

			selectionType = match.arg(testType);
			elimination.pValue <- pvalue; 	# To test if the reduced model is inferior to the full model
			pvalue = pvalue/2;  			# To test that the feature adds statistical power to the model
			
			univariate <- univariateRankVariables(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Kendall",cateGroups,raw.dataFrame,description="Description")
			varMax <- table(univariate$kendall.p<filter.p.value)['TRUE'] + 1;			
			cat("Var Max: ",varMax,"\n");
			if (folds>1)
			{
				cvObject <- crossValidationFeatureSelection(varMax,fraction,pvalue,loops,covariates,Outcome,univariate,data,maxTrainModelSize,
				type,timeOutcome,selectionType,loop.threshold,startOffset,backBootLoops=elimination.bootstrap.steps,
				trainFraction,trainRepetition,elimination.pValue,CVfolds=folds,bootEstimations=bootstrap.steps,interaction,nk=k);
				firstModel <- cvObject$varIDISelection;
				UpdatedModel <- cvObject$updateIDISelection;
				reducedModel <- cvObject$backIDIElimination;
				bootstrappedModel <- cvObject$FullModel.bootstrapped;
			}
			else
			{
				firstModel <- ReclassificationFRESA.Model(varMax,fraction,pvalue,loops,covariates,Outcome,univariate,dataframe=data,maxTrainModelSize,type=type,timeOutcome=timeOutcome,selectionType=selectionType,loop.threshold=loop.threshold,interaction=interaction[1]);
				if (elimination.bootstrap.steps>1)
				{
					cat ("Update\n");
					UpdatedModel <- updateModel(Outcome=Outcome,covariates=covariates,pvalue=c(pvalue,relaxedUpdatePvalue*pvalue),VarFrequencyTable=firstModel$ranked.var,variableListNames=univariate,dataframe=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selectionType,interaction=interaction[2],numberOfModels=1)
					reducedModel <- bootstrapVarElimination(object=UpdatedModel$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=selectionType,bootLoops=elimination.bootstrap.steps,bootFraction=1.0);
				}
				else
				{
					UpdatedModel <- updateModel(Outcome=Outcome,covariates=covariates,pvalue=pvalue,VarFrequencyTable=firstModel$ranked.var,variableListNames=univariate,dataframe=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selectionType,interaction=interaction[2],numberOfModels=1)
					reducedModel <- backVarElimination(object=UpdatedModel$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=selectionType);
				}				
				bootstrappedModel <- bootstrapValidation(fraction,bootstrap.steps,formula(reducedModel$back.formula),Outcome,data,type,plots=FALSE);
			}
		}
		if (OptType == "Residual")
		{	
			elimination.pValue <- pvalue; 	# To test if the reduced model is inferior to the full model
			pvalue = pvalue/2;				# To test that the feature adds statistical power to the model
			testType = match.arg(testType);
			if (testType=="zIDI") 
			{
				testType="Binomial";
			}
			if (length(dependent)==1)
			{
				if (length(table(data[,Outcome]))>2) 
				{
					type = "LM";
					univariate <- univariateRankVariables(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Kendall",cateGroups,raw.dataFrame,description="Description",uniType="Regression")
				}
				else 
				{
					type = "LOGIT";
					univariate <- univariateRankVariables(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Kendall",cateGroups,raw.dataFrame,description="Description")
				}
			}
			else
			{
				univariate <- univariateRankVariables(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Kendall",cateGroups,raw.dataFrame,description="Description")
			}

			bootstrappedModel = NULL;
			varMax <- table(univariate$kendall.p<filter.p.value)['TRUE'] + 1;	
			cat("Var Max: ",varMax,"FitType: ",type," Test Type: ",testType,"\n");
			if (folds>1)
			{
				cvObject <- crossValidationNeRIFeatureSelection(size=varMax,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=univariate,
				dataframe=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,
				loop.threshold=loop.threshold,startOffset=startOffset,trainFraction=trainFraction,trainRepetition=trainRepetition,
				elimination.pValue=elimination.pValue,setIntersect=setIntersect,interaction=interaction,update.pvalue= c(pvalue,relaxedUpdatePvalue*pvalue),backBootLoops=elimination.bootstrap.steps);
				firstModel <- cvObject$FullBaggingModel;
				UpdatedModel <- cvObject$FullUpdatedModel;
				reducedModel <- cvObject$FullBackModel;
			}
			else
			{
				firstModel <- NeRIBasedFRESA.Model(size=varMax,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=univariate,dataframe=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,loop.threshold=loop.threshold,interaction=interaction[1]);
				if (elimination.bootstrap.steps>1)
				{
					UpdatedModel <- updateNeRImodel(Outcome=Outcome,covariates=covariates,pvalue=c(pvalue,relaxedUpdatePvalue*pvalue),VarFrequencyTable=firstModel$ranked.var,variableListNames=univariate,dataframe=data,type=type,testType=testType,timeOutcome=timeOutcome,interaction=interaction[2])
					reducedModel <- bootVarNeRIElimination(object=UpdatedModel$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,bootLoops=elimination.bootstrap.steps,setIntersect=setIntersect);
				}
				else
				{
					UpdatedModel <- updateNeRImodel(Outcome=Outcome,covariates=covariates,pvalue=pvalue,VarFrequencyTable=firstModel$ranked.var,variableListNames=univariate,dataframe=data,type=type,testType=testType,timeOutcome=timeOutcome,interaction=interaction[2])
					reducedModel <- backVarNeRIElimination(object=UpdatedModel$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,setIntersect=setIntersect);
				}
			}
		}
	}
	else
	{
		cat("Expecting a formula object\n");
	}
	result <- list(model = reducedModel$back.model,
	reducedModel = reducedModel,
	univariateAnalysis=univariate,
	firstModel=firstModel,
	updatedModel=UpdatedModel,
	bootstrappedModel=bootstrappedModel,
	cvObject=cvObject,
	used.variables=varMax);
	return (result);
}

