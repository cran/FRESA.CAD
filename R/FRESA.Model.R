FRESA.Model <-
function(formula,data,OptType=c("Binary","Residual"),pvalue=0.05,filter.p.value=0.10,loops=1,maxTrainModelSize=10,loop.threshold=20,elimination.bootstrap.steps=100,bootstrap.steps=100,interaction=c(1,1),print=TRUE,plots=TRUE,CVfolds=10,repeats=1,nk=0,categorizationType=c("Raw","Categorical","ZCategorical","RawZCategorical","RawTail","RawZTail"),cateGroups=c(0.1,0.9),raw.dataFrame=NULL,var.description=NULL,testType=c("zIDI","zNRI","Binomial","Wilcox","tStudent","Ftest"))
{

	categorizationType <- match.arg(categorizationType);

	
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
	
		
		trainFraction <- 1.0-1.0/CVfolds;
		trainRepetition <- repeats*CVfolds;
		fraction = 1.0;
		varMax = nrow(variables);
		baseModel <- paste(dependentout,"~",covariates);
		cvObject = NULL;
		reducedModel = NULL;
		bootstrappedModel = NULL;
		UpdatedModel = NULL;
		filter.z.value <- abs(qnorm(filter.p.value/2))
		if (OptType == "Binary")
		{

			if (length(dependent)==1)
			{
				type = "LOGIT";
			}

			selectionType = match.arg(testType);
			elimination.pValue <- pvalue; 	# To test if the reduced model is inferior to the full model
			pvalue = pvalue/2;  			# To test that the feature adds statistical power to the model
			if (filter.p.value<0.5)
			{
				unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description")
				univariate <- unirank$orderframe;
				varMax <- table(univariate$ZGLM>filter.z.value)['TRUE'] + 1;			
			}
			else
			{
				unirank <- NULL;
				univariate <- variables;
				varMax <- nrow(variables);
			} 
			cat("Z: ",filter.z.value," Var Max: ",varMax,"\n");
			if (CVfolds>1)
			{
				cvObject <- crossValidationFeatureSelection(varMax,fraction,pvalue,loops,covariates,Outcome,timeOutcome,univariate,data,maxTrainModelSize,
				type,selectionType,loop.threshold,startOffset,elimination.bootstrap.steps,
				trainFraction,trainRepetition,elimination.pValue,CVfolds,bootstrap.steps,interaction,nk,unirank,print=print,plots=plots);
				firstModel <- cvObject$varIDISelection;
				UpdatedModel <- cvObject$updateIDISelection;
				reducedModel <- cvObject$backIDIElimination;
				bootstrappedModel <- cvObject$FullModel.bootstrapped;
			}
			else
			{
				firstModel <- ReclassificationFRESA.Model(varMax,fraction,pvalue,loops,covariates,Outcome,univariate,data,maxTrainModelSize,type,timeOutcome,selectionType,loop.threshold,interaction[1]);
				if (length(firstModel$var.names)>0)
				{
					if (elimination.bootstrap.steps>1)
					{
#						cat ("Update\n");
						UpdatedModel <- updateModel(Outcome=Outcome,covariates=covariates,pvalue=c(pvalue,relaxedUpdatePvalue*pvalue),VarFrequencyTable=firstModel$ranked.var,variableList=univariate,data=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selectionType,interaction=interaction[2],numberOfModels=1)
#						cat ("Elimination \n");
						reducedModel <- bootstrapVarElimination(object=UpdatedModel$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=selectionType,loops=elimination.bootstrap.steps,fraction=1.0,print=print,plots=plots);
						bootstrappedModel <- reducedModel$bootCV;
					}
					else
					{
						UpdatedModel <- updateModel(Outcome=Outcome,covariates=covariates,pvalue=pvalue,VarFrequencyTable=firstModel$ranked.var,variableList=univariate,data=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=selectionType,interaction=interaction[2],numberOfModels=1)
						reducedModel <- backVarElimination(object=UpdatedModel$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=selectionType);
						bootstrappedModel <- bootstrapValidation(fraction,bootstrap.steps,reducedModel$back.formula,Outcome,data,type,plots=plots);
					}				
				}
			}
			if (!is.null(bootstrappedModel) && (print)) s <- summary(bootstrappedModel);
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
			if (filter.p.value<0.5)
			{
				if (length(dependent)==1)
				{
					if (length(table(data[,Outcome]))>2) 
					{
						type = "LM";
						unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description",uniType="Regression")
					}
					else 
					{
						type = "LOGIT";
						unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description")
					}
				}
				else
				{
					unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description")
				}
				univariate <- unirank$orderframe;
				varMax <- table(univariate$ZGLM>filter.z.value)['TRUE'] + 1;	
			}
			else
			{
				if (length(table(data[,Outcome]))>2) 
				{
					type = "LM";
				}
				else 
				{
					type = "LOGIT";
				}
				unirank <- NULL;
				univariate <- variables;
				varMax <- nrow(variables);
			}
			bootstrappedModel = NULL;
			cat("Z: ",filter.z.value," Var Max: ",varMax,"FitType: ",type," Test Type: ",testType,"\n");
			if (CVfolds>1)
			{
				cvObject <- crossValidationNeRIFeatureSelection(size=varMax,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,timeOutcome=timeOutcome,variableList=univariate,
				data=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,
				loop.threshold=loop.threshold,startOffset=startOffset,elimination.bootstrap.steps=elimination.bootstrap.steps,trainFraction=trainFraction,trainRepetition=trainRepetition,
				elimination.pValue=elimination.pValue,setIntersect=setIntersect,interaction=interaction,update.pvalue= c(pvalue,relaxedUpdatePvalue*pvalue),unirank=unirank,print=print,plots=plots);
				firstModel <- cvObject$varNeRISelection;
				UpdatedModel <- cvObject$updateNeRISelection;
				reducedModel <- cvObject$backNeRIElimination;
			}
			else
			{
				firstModel <- NeRIBasedFRESA.Model(size=varMax,fraction=fraction,pvalue=pvalue,loops=loops,covariates=covariates,Outcome=Outcome,variableList=univariate,data=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome,loop.threshold=loop.threshold,interaction=interaction[1]);
				if (length(firstModel$var.names)>0)
				{
					if (elimination.bootstrap.steps>1)
					{
						cat ("Update\n");
						UpdatedModel <- updateNeRIModel(Outcome=Outcome,covariates=covariates,pvalue=c(pvalue,relaxedUpdatePvalue*pvalue),VarFrequencyTable=firstModel$ranked.var,variableList=univariate,data=data,type=type,testType=testType,timeOutcome=timeOutcome,interaction=interaction[2])
						cat ("Elimination\n");
						reducedModel <- bootstrapVarNeRIElimination(object=UpdatedModel$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,loops=elimination.bootstrap.steps,setIntersect=setIntersect,print=print,plots=plots);
						bootstrappedModel <- reducedModel$bootCV;
					}
					else
					{
						UpdatedModel <- updateNeRIModel(Outcome=Outcome,covariates=covariates,pvalue=pvalue,VarFrequencyTable=firstModel$ranked.var,variableList=univariate,data=data,type=type,testType=testType,timeOutcome=timeOutcome,interaction=interaction[2])
						reducedModel <- backVarNeRIElimination(object=UpdatedModel$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,setIntersect=setIntersect);
						bootstrappedModel  <- bootstrapValidationNeRI(fraction,bootstrap.steps,reducedModel$back.formula,Outcome,data,type,plots=plots)
					}
				}
			}
		}
	}
	else
	{
		cat("Expecting a formula object\n");
	}
	if (is.null(reducedModel))
	{
		result <- list(final.model = NULL,
		reducedModel = reducedModel,
		univariateAnalysis=univariate,
		firstModel=firstModel,
		updatedModel=UpdatedModel,
		bootstrappedModel=bootstrappedModel,
		cvObject=cvObject,
		used.variables=varMax);
	}
	else
	{
		result <- list(final.model = reducedModel$back.model,
		reducedModel = reducedModel,
		univariateAnalysis=univariate,
		firstModel=firstModel,
		updatedModel=UpdatedModel,
		bootstrappedModel=bootstrappedModel,
		cvObject=cvObject,
		used.variables=varMax);
	}
	return (result);
}


