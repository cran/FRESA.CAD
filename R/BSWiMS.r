BSWiMS.model <-
function(formula=formula,data=NULL,type=c("Auto","LM","LOGIT","COX"),testType=c("Auto","zIDI","zNRI","Binomial","Wilcox","tStudent","Ftest"),pvalue=0.05,elimination.pValue=0.05,update.pvalue=c(0.05,0.05),variableList=NULL,size=0,loops=32,elimination.bootstrap.steps=200,unitPvalues=NULL,adjsize=0.0,fraction=1.0,maxTrainModelSize=20,maxCycles=10,print=FALSE,plots=FALSE)
{
#	print(colnames(data))
#	print(variableList[,1])
#	cat(size,":",pvalue,":",update.pvalue[1],":",elimination.pValue,"\n")
    fforward.model <- NULL
	fupdate.model <- NULL;
	forward.selection.list <- NULL;
	if (class(formula)=="character")
	{
		baseformula <- formula;
		formula <- formula(formula);
	}
	else
	{
		baseformula <- as.character(formula);
		baseformula <- paste(baseformula[2],"~",baseformula[3]);
	}
	varsmod <- all.vars(formula);
		
	varlist <- attr(terms(formula),"variables")
	termslist <- attr(terms(formula),"term.labels");
	setIntersect <- attr(terms(formula),"intercept");
	if (setIntersect == 0)
	{
		covariates = "0";
	}
	else
	{
		covariates = "1";
	}
	startOffset = length(termslist);
	acovariates <- covariates[1];
	if (length(termslist)>0)
	{
		for (i in 1:length(termslist))
		{
			covariates <- paste(covariates,"+",termslist[i]);
			acovariates <- append(acovariates,termslist[i]);
		}
	}
		
	dependent <- as.character(varlist[[2]])
	unitype ="Regression";
	rankingTest="Ztest"
	timeOutcome = NULL;
	Outcome = NULL;
	if (length(dependent)==3)
	{
		timeOutcome = dependent[2];
		Outcome = dependent[3];
		type = "COX";
		if (testType[1]=="Auto") testType="zIDI";
		unitype ="Binary";
		rankingTest="zIDI";
	}
	else
	{
		Outcome = dependent[1];
		if (type[1] == "Auto")
		{
			if (length(table(data[,Outcome]))>2)
			{
				type = "LM";
				if (testType[1]=="Auto") testType="Ftest";
			}
			else
			{
				if (min(data[,Outcome]) != 0)
				{
					type = "LM";
					if (testType[1]=="Auto") testType="Ftest";
				}
				else
				{
					type = "LOGIT";
					if (testType[1]=="Auto") testType="zIDI";
					unitype ="Binary";
					rankingTest="zIDI";
				}
			}
		}
	}

	unirank <- NULL;
	if (is.null(variableList))
	{
		vnames <- names(data);
		names(vnames) <- names(data);
		namesinformula <- vnames %in% all.vars(formula);
		vnames <- vnames[!namesinformula];
		variableList <- cbind(vnames,vnames)
		colnames(variableList) <- c("Name","Description");
		unirank <- uniRankVar(variableList,formula=baseformula,Outcome=Outcome,data=data,categorizationType="Raw",type=type,rankingTest=rankingTest,uniType=unitype,FullAnalysis=FALSE,acovariates=acovariates,timeOutcome=timeOutcome)
		variableList <- unirank$orderframe;
		unirank <- unirank$orderframe;
		if (is.null(unitPvalues))
		{
			unitPvalues <- (1.0-pnorm(variableList$ZUni));
			names(unitPvalues) <- variableList$Name;
			if (size==0)
			{
				size <- sum(p.adjust(unitPvalues,"BH") <= 4*pvalue) + 1; #four times the pvalue for filtering
				if (size>nrow(variableList)) size <- nrow(variableList);
				if (print) cat(nrow(variableList),": Number of variables to test:",size,"\n");
			}
		}
	}
	invariableList <- variableList;
	if (size==0)
	{
		size <- nrow(variableList);
	}
	if (adjsize==0)
	{
		adjsize <- max(size,ncol(data)-2);
	}
	cat(nrow(variableList),": Number of variables to test:",size,"\n");
	isInferiorCnt <-  0;
	isInferior <- 0;
	cycles <- 0;
	firstModel <- NULL;
	BSWiMS.model <- NULL;
	forward.model <- NULL;
	update.model <- NULL;
	formula.list <- character();
	forward.selection.list <- character();
	adjustPvalues <- NULL;
	if (!is.null(unitPvalues))
	{
		adjustPvalues <- p.adjust(unitPvalues,"BH")
		names(adjustPvalues) <-names(unitPvalues);
	}
	metric <- 0;
	while ((isInferiorCnt<2) && (cycles<maxCycles) && (size>1))
	{
		if ((testType=="zIDI") || (testType=="zNRI") )
		{
			forward.model <- ForwardSelection.Model.Bin(size=size,fraction=fraction,pvalue,loops,acovariates,Outcome,variableList,data,maxTrainModelSize,type,timeOutcome,selectionType=testType);
			update.model <- updateModel.Bin(Outcome=Outcome,covariates=covariates,pvalue=update.pvalue,VarFrequencyTable=forward.model$ranked.var,variableList=variableList,data=data,type=type,lastTopVariable= 0,timeOutcome=timeOutcome,selectionType=testType)
			BSWiMS.model <- bootstrapVarElimination_Bin(object=update.model$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=testType,loops=elimination.bootstrap.steps,fraction=fraction,print=print,plots=plots,adjsize=adjsize,uniAdjPvalues=adjustPvalues);
			if (length(forward.model$var.names)>0)
			{
				if (elimination.bootstrap.steps>1)
				{
					currentMeanAUC <- (BSWiMS.model$bootCV$blind.sensitivity + BSWiMS.model$bootCV$blind.specificity)/2;
					metric <- currentMeanAUC
					if (!is.null(firstModel))
					{
						if (length(attr(terms(formula(BSWiMS.model$back.formula)),"term.labels"))>0)
						{
#							firstMeanAUC <- (firstModel$bootCV$blind.sensitivity + firstModel$bootCV$blind.specificity)/2;
							firstAUC <-  (firstModel$bootCV$sensitivity + firstModel$bootCV$specificity)/2;
							curAUC <-  (BSWiMS.model$bootCV$sensitivity + BSWiMS.model$bootCV$specificity)/2;
							firstMedAUC <- 0.85*median(firstAUC);
							if (firstMedAUC < 0.5) firstMedAUC <- 0.5;
							currentMedAUC <- 1.15*median(curAUC);
							if (currentMedAUC>1) currentMedAUC <- 1;
							firstCount <- sum(currentMedAUC >= firstAUC)
							curCount <- sum(firstMedAUC <= curAUC)
							infraction <- 1.0-0.5*(firstCount+curCount)/elimination.bootstrap.steps;
							if (print) 
							{
								cat(BSWiMS.model$back.formula,": Base: ",firstMedAUC,"Current: ",metric," Inferior Count:",firstCount," Tests:",length(firstAUC)," Fraction:",infraction,"\n");
							}
							isInferior <- 1*(infraction>0.90)+1*(infraction>0.999);
						}
						else
						{
							isInferior <-  2;
						}
					}
				}
				else
				{
					BSWiMS.model <- backVarElimination_Bin(object=update.model$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,selectionType=testType,adjsize=adjsize);
					isInferior <-  2;
				}				
			}
		}
		else
		{
#			cat("Forwad\n")
			forward.model <- ForwardSelection.Model.Res(size=size,fraction=fraction,pvalue=pvalue,loops=loops,covariates=acovariates,Outcome=Outcome,variableList=variableList,data=data,maxTrainModelSize=maxTrainModelSize,type=type,testType=testType,timeOutcome=timeOutcome);
#			cat("Update\n")
			update.model <- updateModel.Res(Outcome=Outcome,covariates=covariates,pvalue=update.pvalue,VarFrequencyTable=forward.model$ranked.var,variableList=variableList,data=data,type=type,testType=testType,timeOutcome=timeOutcome)
#			cat("Elimination\n")
			BSWiMS.model <- bootstrapVarElimination_Res(object=update.model$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,loops=elimination.bootstrap.steps,fraction=fraction,setIntersect=setIntersect,print=print,plots=plots,adjsize=adjsize,uniAdjPvalues=adjustPvalues);
			if (length(forward.model$var.names)>0)
			{
				if (elimination.bootstrap.steps>1)
				{
					metric <- BSWiMS.model$bootCV$testRMSE
					if (!is.null(firstModel))
					{
						if (length(attr(terms(formula(BSWiMS.model$back.formula)),"term.labels"))>0)
						{
							firstMedRMS <- 1.15*median(firstModel$bootCV$testSampledRMSE);
							curMedRMS <- 0.85*median(BSWiMS.model$bootCV$testSampledRMSE);
							firstCount <- sum(curMedRMS <= firstModel$bootCV$testSampledRMSE);
							curCount <- sum(firstMedRMS >= BSWiMS.model$bootCV$testSampledRMSE);
							infraction <- 1.0-0.5*(firstCount+curCount)/elimination.bootstrap.steps;
							if (print) 
							{
								cat(BSWiMS.model$back.formula,": Base: ",firstModel$bootCV$testRMSE,"(",max(firstModel$bootCV$testSampledRMSE),") Current: ",BSWiMS.model$bootCV$testRMSE,"(",min(BSWiMS.model$bootCV$testSampledRMSE),") Inferior Count:",firstCount," Tests:",length(firstModel$bootCV$testSampledRMSE)," Fraction:",infraction,"\n");
							}
							isInferior <- 1*(infraction>0.90)+1*(infraction>0.999);
						}
						else
						{
							isInferior <-  2;
						}
					}
				}
				else
				{
					BSWiMS.model <- backVarElimination_Res(object=update.model$final.model,pvalue=elimination.pValue,Outcome=Outcome,data=data,startOffset=startOffset,type=type,testType=testType,setIntersect=setIntersect,adjsize=adjsize);
					isInferior <-  2;
				}
			}
		}
		if (isInferior==0) #removing the models variables
		{
			cat(cycles,":",size,":",metric,":",BSWiMS.model$back.formula,"\n");
#			print(colnames(data))
			formula.list <- append(formula.list,BSWiMS.model$back.formula);
			forward.selection.list <- append(forward.selection.list,forward.model$formula.list);
			termslist <- attr(terms(formula(BSWiMS.model$back.formula)),"term.labels");
			included <- rownames(variableList) %in% termslist;
			variableList <- variableList[!included,]
			size <- size - length(termslist);
			size <- min(size,nrow(variableList));
			isInferiorCnt <- 0;
		}
		if (is.null(firstModel))
		{
			firstModel <- BSWiMS.model;
			fforward.model <- forward.model;
			fupdate.model <- update.model;
			isInferior <- 2*(length(attr(terms(formula(BSWiMS.model$back.formula)),"term.labels"))==0)
		}
		if (isInferior>0)
		{
			isInferiorCnt <- isInferiorCnt+isInferior;
		}
		cycles <- cycles +1;
		isInferior <-  2;
	}
#	print(formula.list);
	if(is.null(unirank))
	{
#		print(invariableList[1:20,])
		unirank <- invariableList;
		unirank$ZUni <- (nrow(invariableList):1)
#		print(unirank[1:20,])
	}
	bagg <- baggedModel(formula.list,data,type,univariate=unirank,useFreq=FALSE); 

	result <- list(BSWiMS.model=firstModel,
		forward.model=fforward.model,
		update.model=fupdate.model,
		univariate=unirank,
		bagging=bagg,
		formula.list=formula.list,
		forward.selection.list=forward.selection.list
	);
	
	return (result);
}