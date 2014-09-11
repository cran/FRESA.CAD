bootVarNeRIElimination <- function (object,pvalue=0.05,Outcome="Class",data,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),bootLoops=250,bootFraction=1.00,setIntersect=1) 
{
  	testType <- match.arg(testType)

	boot.var.NeRISelection <- function (object,pvalue=0.05,Outcome="Class",dataframe,startOffset=0, type = c("LOGIT", "LM","COX"),testType=c("Binomial","Wilcox","tStudent","Ftest"),bootLoops,bootFraction=1.0,setIntersect=1) 
	{
		testType <- match.arg(testType)
		type <- match.arg(type);
		fullModel=object;
		varsList <- as.list(attr(terms(object),"variables"))
		
		climpvalue = pvalue;
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
		NeRICV <- bootstrapValidationNeRI(bootFraction,bootLoops,ftmp,Outcome,dataframe,type,plots=FALSE)
		wcoef <- object$coefficients;
		startSearch = startlist + startOffset;
		if ((length(varsList)-startSearch)>1)
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
				{ # reduce probability by two to test for equivalence of reduced model to the full model
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
					cmax = max(ci[idx],ci2[idx]);
					wcoef[idlist] =  ci2[3];
					if (cmax >= maxPvalue)
					{
						maxPvalue = cmax;
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
				}
			}
	#		cat ("Formula: ",frm1,"\n")
			ftmp <- formula(frm1);
			fullModel <- modelFitting(ftmp,dataframe,type)
			backfrm <- frm1
			if (inherits(fullModel, "try-error"))
			{
				cat("Error: Reduced Formula: ",frm1,"\n");
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
				fullModel <- modelFitting(ftmp,dataframe,type)
			}

		}
		result <- list(Model=fullModel,Removed=removeID,BootModel=NeRICV,backfrm=backfrm);

		return (result)
	}


	changes=1;
	loops=0;
    model <- object;
	while ((changes>0) && (loops<100))
	{
		bk <- boot.var.NeRISelection(object=model,pvalue=pvalue,Outcome=Outcome,dataframe=data,
		startOffset=startOffset,type=type,testType=testType,bootLoops=bootLoops,bootFraction=bootFraction,setIntersect=setIntersect);
		if (!inherits(bk$Model, "try-error"))
		{
			changes = as.integer(bk$Removed);
			if (changes>0)
			{
			  loops = loops + 1;      
			}
			model = bk$Model;
		}
		else
		{
			changes = 1;
			loops = loops + 1;      
		}
	}
	modelReclas <- getVarNeRI(model,dataframe=data,Outcome=Outcome,type);
	if ((bootFraction<1) ||  (changes <0 ))
	{
		NeRICV <- bootstrapValidationNeRI(1.0,bootLoops,model$formula,Outcome,data,type,plots=FALSE);
	}
	else
	{
		NeRICV <- bk$BootModel;
	}
	cat ("Reduced Model\n")
	print(summary(NeRICV$boot.model));

	result <- list(back.model=NeRICV$boot.model,
	loops=loops,
	reclas.info=modelReclas,
	bootCV=NeRICV,
	back.formula=bk$backfrm,
	lastRemoved=changes);
	return (result);
}