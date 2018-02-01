FRESA.Model <-
function(formula,data,OptType=c("Binary","Residual"),pvalue=0.05,filter.p.value=0.10,loops=32,maxTrainModelSize=20,elimination.bootstrap.steps=100,bootstrap.steps=100,print=FALSE,plots=FALSE,CVfolds=1,repeats=1,nk=0,categorizationType=c("Raw","Categorical","ZCategorical","RawZCategorical","RawTail","RawZTail","Tail","RawRaw"),cateGroups=c(0.1,0.9),raw.dataFrame=NULL,var.description=NULL,testType=c("zIDI","zNRI","Binomial","Wilcox","tStudent","Ftest"),lambda="lambda.1se",equivalent=FALSE,bswimsCycles=10,usrFitFun=NULL)
{

	categorizationType <- match.arg(categorizationType);
	cl <- match.call();
	
	cvObject <-  NULL;
	univariate <- NULL;
	eq=NULL;
	bagg=NULL;

	type = "LM";
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
		
		type = "LM";
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
#		cat (termslist,"\n")
		acovariates <- covariates[1];
		if (length(termslist)>0)
		{
			for (i in 1:length(termslist))
			{
				covariates <- paste(covariates,"+",termslist[i]);
				acovariates <- append(acovariates,termslist[i]);
			}
		}
#		cat (covariates,"\n")
#		print(acovariates);

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
	
		if (CVfolds>nrow(data))
		{
			cat("Setting to LOO CV\n");
			CVfolds=nrow(data);
		}
		
		trainFraction <- 1.0-1.0/CVfolds;
		trainRepetition <- repeats*CVfolds;
		fraction = 1.0000;					# will be working with 1.0000 fraction of the samples for bootstrap training 
		varMax = nrow(variables);
		baseModel <- paste(dependentout,"~",covariates);
		cvObject = NULL;
		reducedModel = NULL;
		bootstrappedModel = NULL;
		UpdatedModel = NULL;
		filter.z.value <- abs(qnorm(filter.p.value/2))
		cutpvalue <- 4.0*filter.p.value
		if (cutpvalue>0.5) cutpvalue=0.5;
		selectionType = match.arg(testType);
		testType = match.arg(testType);
		if (((length(table(data[,Outcome]))>2)||(min(data[,Outcome])<0))&&(OptType == "Binary"))
		{
			OptType = "Residual";
		}
		adjsize <- 1.0;
		if (OptType == "Binary")
		{


			if (length(dependent)==1)
			{
				type = "LOGIT";
			}

			elimination.pValue <- pvalue; 	# To test if the variable is part of the model 
			unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="zIDI",cateGroups,raw.dataFrame,description="Description",uniType="Binary",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome);
			univariate <- unirank$orderframe;
			unitPvalues <- (1.0-pnorm(univariate$ZUni));
			names(unitPvalues) <- univariate$Name;
			
#			print(univariate[1:20,]);
			varMax <- sum(univariate$ZUni >= filter.z.value) + 1;
			pvarMax <- sum(p.adjust(unitPvalues,"BH") < 2*filter.p.value)+1;
			unitPvalues[unitPvalues>0.5] <- 0.5;
			cat("Unadjusted size:",varMax," Adjusted Size:",pvarMax,"\n")
			varMax <- min(c(pvarMax,varMax));
			if (is.null(varMax)) varMax <- 2;
			if (is.na(varMax)) varMax <- 2;
			if (is.nan(varMax)) varMax <- 2;
			if (varMax >  nrow(variables)) varMax = nrow(variables);

			redlist <- 2*(1.0-pnorm(univariate$ZUni)) < cutpvalue;
			totlist <- min(sum(1*redlist),100);
			
			if (categorizationType == "RawRaw")  
			{
				sdata <- data;
				reddata <- data;
			}
			else
			{
				if (is.na(timeOutcome))
				{
					reddata <- data[,c(Outcome,acovariates[-1],unique(as.character(univariate[redlist,2])))];
					sdata <- data[sample(nrow(data),min(nrow(data),100)),c(Outcome,acovariates[-1],unique(as.character(univariate[1:totlist,2])))];
				}
				else
				{
					reddata <- data[,c(Outcome,timeOutcome,acovariates[-1],unique(as.character(univariate[redlist,2])))];
					sdata <- data[sample(nrow(data),min(nrow(data),100)),c(Outcome,timeOutcome,acovariates[-1],unique(as.character(univariate[1:totlist,2])))];
				}
			}

			randomModel <- ForwardSelection.Model.Res(totlist,1.0,0.05,250,"1",Outcome,univariate,sdata,maxTrainModelSize,"LM","Ftest",timeOutcome,1,randsize= -1);
			adjsize <- as.integer(min(40*randomModel$random.formula.size,ncol(data)*ncol(data)));
			if (adjsize<1) adjsize=1;
		
			
			cat("\n Z: ",filter.z.value,", Var Max: ",varMax,", s1:",ncol(data),", s2:",ncol(reddata),", Independent Size:",adjsize,"\n");
			shorUniv <- univariate[redlist,]

			
			if (CVfolds>1)
			{
				if (categorizationType!="RawRaw")  
				{
					rownames(variables) <- variables[,1];
					unirank$variableList <- variables[unique(as.character(univariate[redlist,2])),]
				}
#				BSWiMS.models <- BSWiMS.model(formula=formula,data=reddata,type=type,testType=selectionType,pvalue=pvalue,elimination.pValue=elimination.pValue,update.pvalue=c(pvalue,pvalue),variableList=shorUniv,size=varMax,loops=loops,elimination.bootstrap.steps=bootstrap.steps,unitPvalues=unitPvalues,adjsize=adjsize,fraction=1.0,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots);
				cvObject <- crossValidationFeatureSelection_Bin(varMax,fraction,pvalue,loops,acovariates,Outcome,timeOutcome,NULL,reddata,maxTrainModelSize,type,selectionType,startOffset,elimination.bootstrap.steps,trainFraction,trainRepetition,elimination.pValue,bootstrap.steps,nk,unirank,print=print,plots=plots,lambda=lambda,adjsize=adjsize,equivalent=equivalent,bswimsCycles=bswimsCycles,usrFitFun);
				firstModel <- cvObject$forwardSelection;
				UpdatedModel <- cvObject$updateforwardSelection;
				reducedModel <- cvObject$BSWiMS;
				bootstrappedModel <- cvObject$FullBSWiMS.bootstrapped;
				BSWiMS.models <- cvObject$BSWiMS.models;
			}
			else
			{
				BSWiMS.models <- BSWiMS.model(formula=formula,data=reddata,type=type,testType=selectionType,pvalue=pvalue,elimination.pValue=elimination.pValue,update.pvalue=c(pvalue,pvalue),variableList=shorUniv,size=varMax,loops=loops,elimination.bootstrap.steps=bootstrap.steps,unitPvalues=unitPvalues,adjsize=adjsize,fraction=1.0,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots);
				firstModel <- BSWiMS.models$forward.model;
				UpdatedModel <- BSWiMS.models$update.model;
				reducedModel <- BSWiMS.models$BSWiMS.model;
				bootstrappedModel <- reducedModel$bootCV;
			
			}

		}
		if (OptType == "Residual")
		{	
			elimination.pValue <- pvalue; 	# To test if the variable is part of the model
			if (testType=="zIDI") 
			{
				testType="Ftest";
			}
			if (length(dependent)==1)
			{
				if ((length(table(data[,Outcome]))>2)||(min(data[,Outcome])<0))
				{
					type = "LM";
					unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description",uniType="Regression",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
				}
				else 
				{
					if (type == "LM") type = "LOGIT";
					unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description",uniType="Binary",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
				}
			}
			else
			{
				    unirank <- uniRankVar(variables,baseModel,Outcome,data,categorizationType,type,rankingTest="Ztest",cateGroups,raw.dataFrame,description="Description",uniType="Binary",FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
			}
			univariate <- unirank$orderframe;
			unitPvalues <- (1.0-pnorm(univariate$ZUni));
			unitPvalues[unitPvalues>0.5] <- 0.5;
			names(unitPvalues) <- univariate$Name;

			varMax <- sum(univariate$ZUni >= filter.z.value) + 1;	
			pvarMax <- sum(p.adjust(unitPvalues,"BH") < 4*filter.p.value)+1;
			cat("Unadjusted size:",varMax," Adjusted Size:",pvarMax,"\n")
			varMax <- min(c(pvarMax,varMax));
			if (is.null(varMax)) varMax <- 2;
			if (is.na(varMax)) varMax <- 2;
			if (is.nan(varMax)) varMax <- 2;
			if (varMax >  nrow(variables)) varMax = nrow(variables);
			bootstrappedModel = NULL;
			redlist <- (2*(1.0-pnorm(univariate$ZUni)) < cutpvalue);


			totlist <- min(sum(1*redlist),100);

			if (categorizationType=="RawRaw")  
			{
				sdata <- data;
				reddata <- data;
			}
			else
			{
				if (is.na(timeOutcome))
				{
					reddata <- data[,c(Outcome,acovariates[-1],unique(as.character(univariate[redlist,2])))];
					sdata <- data[sample(nrow(data),min(nrow(data),100)),c(Outcome,acovariates[-1],unique(as.character(univariate[1:totlist,2])))];
				}
				else
				{
					reddata <- data[,c(Outcome,timeOutcome,acovariates[-1],unique(as.character(univariate[redlist,2])))];
					sdata <- data[sample(nrow(data),min(nrow(data),100)),c(Outcome,timeOutcome,acovariates[-1],unique(as.character(univariate[1:totlist,2])))];
				}
			}

			randomModel <- ForwardSelection.Model.Res(totlist,1.0,0.05,250,"1",Outcome,univariate,sdata,maxTrainModelSize,"LM","Ftest",timeOutcome,1,randsize= -1);
			adjsize <- as.integer(min(40*randomModel$random.formula.size,ncol(data)*ncol(data)));

			if (adjsize<1) adjsize=1;
			cat("\n Z: ",filter.z.value," Var Max: ",varMax,"FitType: ",type," Test Type: ",testType,"s1:",ncol(data),"s2:",ncol(reddata),", Independent Size:",adjsize,"\n");
			shorUniv <- univariate[redlist,]

			if (CVfolds>1)
			{
				if (categorizationType != "RawRaw")  
				{
					rownames(variables) <- variables[,1];
					unirank$variableList <- variables[unique(as.character(univariate[redlist,2])),]
				}
				cvObject <- crossValidationFeatureSelection_Res(size=varMax,fraction=fraction,pvalue=pvalue,loops=loops,covariates=acovariates,Outcome=Outcome,timeOutcome=timeOutcome,variableList=unirank$variableList,
				data=reddata,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,startOffset=startOffset,elimination.bootstrap.steps=elimination.bootstrap.steps,trainFraction=trainFraction,trainRepetition=trainRepetition,
				elimination.pValue=elimination.pValue,setIntersect=setIntersect,update.pvalue= c(pvalue,pvalue),unirank=unirank,print=print,plots=plots,lambda=lambda,adjsize=adjsize,equivalent=equivalent,bswimsCycles=bswimsCycles,usrFitFun=usrFitFun);
				firstModel <- cvObject$forwardSelection;
				UpdatedModel <- cvObject$updatedforwardModel;
				reducedModel <- cvObject$BSWiMS;
				bootstrappedModel <- cvObject$BSWiMS$bootCV;
				BSWiMS.models <- cvObject$BSWiMS.models;
			}
			else
			{
				BSWiMS.models <- BSWiMS.model(formula=formula,data=reddata,type=type,testType=testType,pvalue=pvalue,elimination.pValue=elimination.pValue,update.pvalue=c(pvalue,pvalue),variableList=shorUniv,size=varMax,loops=loops,elimination.bootstrap.steps=bootstrap.steps,unitPvalues=unitPvalues,adjsize=adjsize,fraction=1.0,maxTrainModelSize=maxTrainModelSize,maxCycles=bswimsCycles,print=print,plots=plots);
				firstModel <- BSWiMS.models$forward.model;
				UpdatedModel <- BSWiMS.models$update.model;
				reducedModel <- BSWiMS.models$BSWiMS.model;
				bootstrappedModel <- reducedModel$bootCV;

			}
		}
	}
	else
	{
		cat("Expecting a formula object\n");
	}
	if (is.null(reducedModel))
	{
		result <- list(BSWiMS.model = NULL,
		reducedModel = reducedModel,
		univariateAnalysis=univariate,
		forwardModel=firstModel,
		updatedforwardModel=UpdatedModel,
		bootstrappedModel=bootstrappedModel,
		cvObject=cvObject,
		used.variables=varMax,
		independenSize=adjsize,
		call=cl);
	}
	else
	{
	
		collectFormulas <- BSWiMS.models$forward.selection.list
		bagg <- baggedModel(collectFormulas,data,type,Outcome,timeOutcome,univariate=univariate,useFreq=loops);
		eq <- NULL;		
		if (length(reducedModel$back.model$coefficients) > 1 ) 
		{
			shortcan <- bagg$frequencyTable[(bagg$frequencyTable >= (loops*0.05))];
			eadjsize <- length(unitPvalues);
			if (adjsize<eadjsize) unitPvalues <- unitPvalues[1:adjsize];
#			cat(adjsize,":",length(unitPvalues),"\n")
			apvalues <- p.adjust(unitPvalues,"bonferroni",adjsize);
			modeltems <- attr(terms(reducedModel$back.model),"term.labels");
			univariatenames <- names(apvalues[apvalues < pvalue/length(modeltems)]);
#			print(univariatenames);
			eshortlist <- unique(c(names(shortcan),str_replace_all(modeltems,":","\\*")));
			eshortlist <- unique(c(eshortlist,univariatenames));
			eshortlist <- eshortlist[!is.na(eshortlist)];
#			print(eshortlist);
			if (length(eshortlist)>0)
			{
#				print(eshortlist);
				nameslist <- c(all.vars(BSWiMS.models$bagging$bagged.model$formula),as.character(univariate[eshortlist,2]));
				nameslist <- unique(nameslist[!is.na(nameslist)]);
#				print(nameslist);
				if (categorizationType != "RawRaw") 
				{
					eqdata <- data[,nameslist];
				}
				else
				{
					eqdata <- data;
				}
				eq <- reportEquivalentVariables(reducedModel$back.model,pvalue = pvalue,
							  data=eqdata,
							  variableList=cbind(eshortlist,eshortlist),
							  Outcome = Outcome,
							  timeOutcome=timeOutcome,								  
							  type = type,
							  eqPGain = 100.00,method="bonferroni",osize=adjsize);
			}
		}

		result <- list(BSWiMS.model = BSWiMS.models$bagging$bagged.model,
		reducedModel = reducedModel,
		univariateAnalysis=univariate,
		forwardModel=firstModel,
		updatedforwardModel=UpdatedModel,
		bootstrappedModel=bootstrappedModel,
		cvObject=cvObject,
		used.variables=varMax,
		independenSize=adjsize,
		bagging=bagg,
		eBSWiMS.model=eq,
		BSWiMS.models=BSWiMS.models,
		call=cl
		);
	}
	return (result);
}


