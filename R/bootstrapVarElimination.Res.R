bootstrapVarElimination_Res <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loops=250,fraction=1.00,setIntersect=1,print=TRUE,plots=TRUE,adjsize=1) 
{
  	testType <- match.arg(testType)

	boot.var.NeRISelection <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),loops,fraction=1.0,setIntersect=1,NeRICV=NULL) 
	{
		testType <- match.arg(testType)
		type <- match.arg(type);
		FullModel=object;
		varsList <- as.list(attr(terms(object),"variables"))
		
		climpvalue = max(5*pvalue,0.35); # for % 
		removeID = 0;

		outCome = paste(varsList[2]," ~ ",setIntersect);
		startlist = 3 ;
		frm1 = outCome;
		for ( i in startlist:length(varsList))
		{
			frm1 <- paste(frm1,paste(" + ",varsList[i]));
		}
		ftmp <- formula(frm1);
		backfrm <- frm1;
#		cat("Start  Formula :",frm1,"\n")
		if (is.null(NeRICV)) NeRICV <- bootstrapValidation_Res(fraction,loops,ftmp,Outcome,data,type,plots=plots)
#		cat("Bootstrapped  Formula :",frm1,"\n")
		startSearch = startlist + startOffset;
		maxPvalue=pvalue;
		if (length(varsList)>startSearch)
		{		
			idlist=startOffset+1;
			frm1 = outCome;
			if ((startSearch-1) >= startlist)
			{
				for ( i in startlist:(startSearch-1))
				{
					frm1 <- paste(frm1,paste(" + ",varsList[i]));
				}
			}
	#		cat("Elimination Base Formula :",frm1,"\n")
			idx = 2;
			who = -1;
			maxPvalue = pvalue;
			for ( i in startSearch:length(varsList))
			{
				if (any(is.na(NeRICV$bin.pvlaues[,idlist])))
				{
					who = i;
				}
				else
				{ # reduce probability by two to test for equivalence of reduced model to the Full model
					switch(testType, 
						tStudent = 
						{ 
							ci <- as.vector(quantile(NeRICV$tStudent.pvalues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
							ci2 <- as.vector(quantile(NeRICV$test.tStudent.pvalues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
						},
						Wilcox = 
						{ 
							ci <- as.vector(quantile(NeRICV$wilcox.pvalues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
							ci2 <- as.vector(quantile(NeRICV$test.wilcox.pvalues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
						},
						Binomial =
						{ 
							ci <- as.vector(quantile(NeRICV$bin.pvlaues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
							ci2 <- as.vector(quantile(NeRICV$test.bin.pvlaues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
						},
						Ftest =
						{ 
							ci <- as.vector(quantile(NeRICV$F.pvlaues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
							ci2 <- as.vector(quantile(NeRICV$test.F.pvlaues[,idlist], probs = c(climpvalue, 0.5, 1-climpvalue), na.rm = TRUE,names = FALSE, type = 7));
						},
					)
					cmax = max(ci[2],ci2[2]);
					if (cmax >= maxPvalue)
					{
						maxPvalue = cmax;
						who = i;
					}
					if ((ci[3] > pvalue)&&(maxPvalue == pvalue))
					{
						who = i;
					}
				}
				idlist=idlist+1;
			}
			if ((length(varsList) == startSearch) && (who == startSearch)) 
			{
				who = -1;
				removeID = -1;
			}
			for ( i in startSearch:length(varsList))
			{
				if (who != i)
				{
					frm1 <- paste(frm1,paste(" + ",varsList[i]));
				}
				else
				{
					removeID=idlist;
#					cat ("Removed ",paste(" -> ",varsList[i]),"\n");
				}
			}
#			cat ("Rm: Formula: ",frm1,"\n")
			ftmp <- formula(frm1);
			FullModel <- modelFitting(ftmp,data,type)
			backfrm <- frm1
			if (inherits(FullModel, "try-error"))
			{
#				cat("Error: Reduced Formula: ",frm1,"\n");
				who= length(varsList)
				frm1 = outCome;
				if ((startSearch-1) >= startlist)
				{
					for ( i in startlist:(startSearch-1))
					{
						frm1 <- paste(frm1,paste(" + ",varsList[i]));
					}
				}
				for ( i in startSearch:length(varsList))
				{
					if (who != i)
					{
						frm1 <- paste(frm1,paste(" + ",varsList[i]));
					}
				}
				ftmp <- formula(frm1);
				backfrm <- frm1
				FullModel <- modelFitting(ftmp,data,type)
			}
			if (removeID>0)	NeRICV <- bootstrapValidation_Res(fraction,loops,ftmp,Outcome,data,type,plots=plots)


		}
		result <- list(Model=FullModel,Removed=removeID,BootModel=NeRICV,backfrm=backfrm,max.pvalue=maxPvalue);

		return (result)
	}

	bkobj <- NULL;
	minbootRMSE <- 1.0;
	stopFratio = 1.025; # stop if there is a small increase in RMS test error 

	if (adjsize>1)
	{
		bkobj <- bootstrapVarElimination_Res(object,pvalue,Outcome,data,startOffset,type,testType,loops,fraction,setIntersect,print,plots,adjsize=1); #just remove the features that produce similar models
		object <- bkobj$back.model;
		adjsize = floor(adjsize);	
		adjsize <- min(adjsize,ncol(data)-1);
		minbootRMSE <- bkobj$bootCV$testRMSE;
		bkobj <- bkobj$bootCV;
#		cat(" Adjust size:",adjsize,"\n");
		stopFratio = 1.25; # max increase in RMS test error 
	}
	else
	{
		bkobj <- bootstrapValidation_Res(fraction,loops,formula(object),Outcome,data,type,plots=plots)
		minbootRMSE <- bkobj$testRMSE;		
	}
	
	startRMSE <- minbootRMSE;
	varsList <- as.list(attr(terms(object),"variables"))
	outCome = paste(varsList[2]," ~ ",setIntersect);
	startlist = 3 ;
	frm1 = outCome;
	for ( i in startlist:length(varsList))
	{
		frm1 <- paste(frm1,paste(" + ",varsList[i]));
	}
	model.formula <- frm1;
	
	
	changes=1;
	loopsAux=0;
    model <- object;
	beforeFSCmodel <- object;
	beforeFSC.model.formula <- frm1;
	wts <- rep(1,length(beforeFSCmodel$coefficients));
	names(wts) <- names(beforeFSCmodel$coefficients);
	changes2 <- 0

	lastbootval = NULL;
	while ((changes>0) && (loopsAux<100)) 
	{
		p.elimin <- pvalue;
		if (adjsize>1)
		{
			modsize <- length(as.list(attr(terms(model),"term.labels")));
			if (modsize<1) modsize=1;
			p.elimin <- min(pvalue,2*modsize*pvalue/adjsize) # # BH alpha the elimination p-value
		}

		bk <- boot.var.NeRISelection(object=model,pvalue=p.elimin,Outcome=Outcome,data=data,
		startOffset=startOffset,type=type,testType=testType,loops=loops,fraction=fraction,setIntersect=setIntersect,lastbootval);
				
		if (!inherits(bk$Model, "try-error"))
		{
			lastbootval <- bk$BootModel;
			weight <- 1.0;
#			cat ("RMSE :",bk$BootModel$testRMSE,"\n")
			if (bk$BootModel$testRMSE >= minbootRMSE) 
			{
				weight <- 1.0;	# the weight of the old coefficient estimations
			}
			else
			{
				minbootRMSE <- bk$BootModel$testRMSE;
				weight <- 0.0;	# suppress the betas beyond the minimum
			}
			changes = as.integer(bk$Removed);
			if (bk$BootModel$testRMSE > stopFratio*minbootRMSE) 
			{
				changes=0;	# stop with an increase from the minimum test error
			}
			if (changes>0)
			{
				loopsAux = loopsAux + 1;
				changes2<- as.character(as.list(attr(terms(model),"variables")))[which(!(as.character(as.list(attr(terms(model),"variables")))%in%as.character(as.list(attr(terms(bk$Model),"variables")))))]
				model = bk$Model;
				model.formula <- bk$backfrm;
				if (length(changes2)>1)
				{
					changes2<-changes2[2]
				}
				if (adjsize>1)
				{
					for (i in 1:length(beforeFSCmodel$coefficients))
					{
						notadded = TRUE;
						for (j in 1:length(model$coefficients))
						{
							if (names(beforeFSCmodel$coefficients)[i] == names(model$coefficients)[j])
							{
								beforeFSCmodel$coefficients[i] <- (weight*beforeFSCmodel$coefficients[i] + (1.0-weight)*model$coefficients[j]); # it will weight based on the probability
								notadded=FALSE;
							}
						}
						if (notadded)
						{
							beforeFSCmodel$coefficients[i] <- weight*beforeFSCmodel$coefficients[i]; # it will average with zero
							wts[i] = wts[i]*weight;
						}
					}
				}
			}
			if (changes < 0)
			{
				changes2<- changes
			}
		}
		else
		{
			changes = 1;
			loopsAux = loopsAux + 1;      
		}
#		model = bk$Model;
	}
	modelReclas <- getVar.Res(model,data=data,Outcome=Outcome,type);
	NeRICV <- bootstrapValidation_Res(fraction,loops,model.formula,Outcome,data,type,plots=plots);
	if ((adjsize>1)&&(print==TRUE)) 
	{
		cat("Before BSC Mod:",beforeFSC.model.formula,"\n");
		cat("Final  Formula:",model.formula,"\n")
		cat("Adjust size:",adjsize,"\n");
		cat("Start RMSE:",startRMSE,"final RSME:",NeRICV$testRMSE,"Min RMSE:",minbootRMSE,"\n")
#		print(wts,digits=3)
#		cat("\n");
	}
#	print(summary(NeRICV$boot.model));

	result <- list(back.model=model,
	loops=loopsAux,
	reclas.info=modelReclas,
	bootCV=NeRICV,
	back.formula=formula(model.formula),
	lastRemoved=changes2,
	beforeFSC.model = beforeFSCmodel,
	beforeFSC.formula = formula(beforeFSC.model.formula));
	return (result);
}